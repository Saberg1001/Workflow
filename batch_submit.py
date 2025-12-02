import argparse
import csv
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from smiles_to_input import repalce_coordinates, smiles_to_gjf

JOB_LOG_DEFAULT = Path("jobs.csv")
TERMINAL_STATES = {
    "COMPLETED",
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
    "OUT_OF_MEMORY",
    "PREEMPTED",
}
JOB_FIELDS = [
    "smile",
    "work_dir",
    "gjf",
    "script",
    "job_id",
    "status",
    "submitted_at",
    "completed_at",
]


def load_smiles(smiles_file: Path) -> List[str]:
    lines = smiles_file.read_text(encoding="utf-8").splitlines()
    smiles = [
        line.split("#", 1)[0].strip()
        for line in lines
        if line.strip() and not line.strip().startswith("#")
    ]
    if not smiles:
        raise ValueError(f"No SMILES found in {smiles_file}")
    return smiles


def prepare_jobs(
    smiles: Iterable[str],
    template: Path,
    submit_template: Path,
) -> List[Dict[str, str]]:
    jobs: List[Dict[str, str]] = []
    for smile in smiles:
        work_dir, gjf_dir = smiles_to_gjf(smile)
        repalce_coordinates(work_dir, gjf_dir, smile, template, submit_template)
        job_script = work_dir / f"{smile}.sh"
        jobs.append(
            {
                "smile": smile,
                "work_dir": str(work_dir.resolve()),
                "gjf": str(gjf_dir.resolve()),
                "script": str(job_script.resolve()),
            }
        )
    return jobs


def append_jobs_to_log(log_path: Path, records: List[Dict[str, str]]) -> None:
    new_file = not log_path.exists()
    with log_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=JOB_FIELDS)
        if new_file:
            writer.writeheader()
        for record in records:
            writer.writerow(record)


def read_job_log(log_path: Path) -> List[Dict[str, str]]:
    if not log_path.exists():
        return []
    with log_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def run_command(args: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        args,
        capture_output=True,
        text=True,
        check=False,
    )


def submit_job(script_path: Path) -> str:
    result = run_command(["sbatch", str(script_path)])
    if result.returncode != 0:
        raise RuntimeError(
            f"sbatch failed for {script_path}: {result.stderr.strip() or result.stdout.strip()}"
        )
    parts = result.stdout.strip().split()
    if not parts:
        raise RuntimeError(f"Unexpected sbatch output: {result.stdout}")
    job_id = parts[-1]
    if not job_id.isdigit():
        raise RuntimeError(f"Could not parse job id from sbatch output: {result.stdout}")
    return job_id


def refresh_job_status(job: Dict[str, str]) -> None:
    job_id = job.get("job_id")
    if not job_id:
        return

    state = ""
    sacct = run_command(
        [
            "sacct",
            "-j",
            job_id,
            "--format=JobID,State",
            "-P",
            "-n",
        ]
    )
    if sacct.returncode == 0 and sacct.stdout.strip():
        for line in sacct.stdout.strip().splitlines():
            fields = line.split("|")
            if not fields:
                continue
            job_identifier = fields[0]
            job_state = fields[1] if len(fields) > 1 else ""
            if job_identifier == job_id or job_identifier.startswith(f"{job_id}."):
                state = job_state
                break

    if not state:
        squeue = run_command(
            [
                "squeue",
                "-j",
                job_id,
                "-h",
                "-o",
                "%T",
            ]
        )
        if squeue.returncode == 0:
            state = squeue.stdout.strip()

    if not state:
        state = job.get("status") or "UNKNOWN"

    job["status"] = state
    if state in TERMINAL_STATES and not job.get("completed_at"):
        job["completed_at"] = datetime.now(timezone.utc).isoformat()


def overwrite_log(log_path: Path, jobs: List[Dict[str, str]]) -> None:
    with log_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=JOB_FIELDS)
        writer.writeheader()
        for job in jobs:
            writer.writerow(job)


def print_status_table(jobs: List[Dict[str, str]]) -> None:
    columns = ["smile", "job_id", "status", "submitted_at", "completed_at"]
    widths = {col: len(col) for col in columns}
    for job in jobs:
        for col in columns:
            widths[col] = max(widths[col], len(job.get(col, "") or ""))
    header = " | ".join(col.ljust(widths[col]) for col in columns)
    print(header)
    print("-" * len(header))
    for job in jobs:
        row = " | ".join((job.get(col, "") or "").ljust(widths[col]) for col in columns)
        print(row)


def watch_jobs(log_path: Path, interval: int) -> None:
    while True:
        jobs = read_job_log(log_path)
        if not jobs:
            print("No jobs in log; nothing to watch.")
            return
        for job in jobs:
            refresh_job_status(job)
        overwrite_log(log_path, jobs)
        print_status_table(jobs)
        if all(job.get("status") in TERMINAL_STATES for job in jobs):
            print("All jobs reached terminal state.")
            return
        time.sleep(interval)


def handle_submit(args: argparse.Namespace) -> None:
    smiles_file = Path(args.smiles_file)
    template = Path(args.template)
    submit_template = Path(args.submit_script)
    log_path = Path(args.job_log)

    smiles = load_smiles(smiles_file)
    jobs = prepare_jobs(smiles, template, submit_template)

    submission_records: List[Dict[str, str]] = []
    for job in jobs:
        job_id = submit_job(Path(job["script"]))
        now = datetime.now(timezone.utc).isoformat()
        print(f"Submitted {job['smile']} -> Job {job_id}")
        submission_records.append(
            {
                **job,
                "job_id": job_id,
                "status": "SUBMITTED",
                "submitted_at": now,
                "completed_at": "",
            }
        )

    append_jobs_to_log(log_path, submission_records)

    if args.watch:
        watch_jobs(log_path, args.poll_interval)


def handle_status(args: argparse.Namespace) -> None:
    log_path = Path(args.job_log)
    if args.watch:
        watch_jobs(log_path, args.poll_interval)
        return
    jobs = read_job_log(log_path)
    if not jobs:
        print(f"No job log entries in {log_path}")
        return
    for job in jobs:
        refresh_job_status(job)
    overwrite_log(log_path, jobs)
    print_status_table(jobs)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch Gaussian submission and status tracking for Slurm."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    submit = subparsers.add_parser("submit", help="Generate inputs and submit jobs.")
    submit.add_argument(
        "--smiles-file",
        required=True,
        help="Path to a text file containing one SMILES string per line.",
    )
    submit.add_argument(
        "--template",
        default="template.gjf",
        help="Gaussian template to use for link sections.",
    )
    submit.add_argument(
        "--submit-script",
        default="submit.sh",
        help="Submit-script template copied per molecule.",
    )
    submit.add_argument(
        "--job-log",
        default=str(JOB_LOG_DEFAULT),
        help="CSV file used to persist job metadata.",
    )
    submit.add_argument(
        "--watch",
        action="store_true",
        help="Watch submitted jobs until they reach terminal states.",
    )
    submit.add_argument(
        "--poll-interval",
        type=int,
        default=60,
        help="Polling interval in seconds when watching job status.",
    )
    submit.set_defaults(func=handle_submit)

    status = subparsers.add_parser("status", help="Refresh and print job status.")
    status.add_argument(
        "--job-log",
        default=str(JOB_LOG_DEFAULT),
        help="CSV file used to persist job metadata.",
    )
    status.add_argument(
        "--watch",
        action="store_true",
        help="Continuously poll job status until completion.",
    )
    status.add_argument(
        "--poll-interval",
        type=int,
        default=60,
        help="Polling interval in seconds when watching job status.",
    )
    status.set_defaults(func=handle_status)

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main(sys.argv[1:])

