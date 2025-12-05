from pathlib import Path
import subprocess

pair_dirs = Path(".").glob("pair*")

for pair in pair_dirs:
    for sub_dir_name in ["amine", "aldehyde"]:
        sub_dir = pair / sub_dir_name
        submit_script = sub_dir / "submit.sh"

        if not submit_script.is_file():
            continue

        print(f"Submitting in {sub_dir}")

        # 1. 修正行尾
        try:
            subprocess.run(
                ["sed", "-i", "s/\r$//", "submit.sh"],
                cwd=str(sub_dir),
                check=True,
            )
            print("fixed submit.sh")
        except subprocess.CalledProcessError:
            print("ERROR: failed to fix line ending")
            continue  # 行尾修不好就不用再提交了

        # 2. 提交作业
        try:
            subprocess.run(
                ["sbatch", "submit.sh"],
                cwd=str(sub_dir),
                check=True,
            )
            print("submit success")
        except subprocess.CalledProcessError:
            print("ERROR: submit failed")