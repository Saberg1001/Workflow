import subprocess
import sys
from datetime import datetime
from pathlib import Path
import config
import csv
import pandas as pd
import re


def extract_last_task(work_dir, output_name="excited.log"):
    output_log_dir = Path(work_dir / "output.log")
    output_dir = Path(work_dir / output_name)
    if not Path(output_log_dir).exists():
        print(f"错误: 文件 '{work_dir}' 不存在！")
        sys.exit(1)

    try:
        with open(output_log_dir, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        positions = []
        for i, line in enumerate(lines):
            if "Initial command:" in line:
                positions.append(i)

        if not positions:
            print("错误: 未找到 'Initial command:'")
            sys.exit(1)

        last_position = positions[-1]
        extracted_content = lines[last_position:]

        with open(output_dir, "w", encoding="utf-8") as f:
            f.writelines(extracted_content)

        print(f"成功提取最后一个任务到 '{output_dir}'")

    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)
    return output_dir


def get_CDFT_data(sub_dir):
    if not Path(sub_dir).exists():
        msg = f"目录 '{sub_dir}' 不存在！"
        print(f"错误: {msg}")
        return msg, None

    neutral_dir = Path(sub_dir / "neutral.fchk")
    anion_dir = Path(sub_dir / "anion.fchk")
    cation_dir = Path(sub_dir / "cation.fchk")
    missing = [
        str(path.name)
        for path in (neutral_dir, anion_dir, cation_dir)
        if not path.exists()
    ]
    if missing:
        msg = f"缺少 {', '.join(missing)}"
        print(f"错误: {sub_dir} {msg}")
        return msg, None

    input_text = f"""\
        {neutral_dir}
        22
        2
        {neutral_dir}
        {anion_dir}
        {cation_dir}
        q
        """
    input_file = Path(sub_dir / "part1_input.txt")
    with open(input_file, "w", encoding="utf-8") as f:
        f.write(input_text)

    CDFT_log = Path(sub_dir / "CDFT.log")
    CDFT_txt = Path(sub_dir / "CDFT.txt")
    run_status = "failed"

    try:
        with open(input_file, "r", encoding="utf-8") as stdin_file, open(
            CDFT_log, "w", encoding="utf-8"
        ) as CDFT_stream:
            subprocess.run(
                ["Multiwfn_noGUI"],
                stdin=stdin_file,
                stderr=subprocess.STDOUT,
                stdout=CDFT_stream,
                cwd=sub_dir,
                text=True,
            )
    except Exception as e:
        msg = f"运行异常: {e}"
        print(f"错误: {msg}")
        return msg, None

    with open(CDFT_log, "r", encoding="utf-8") as f:
        content = f.read()
        if "Done! Data have been outputted to CDFT.txt in current folder!" in content:
            run_status = "success"
    with open(CDFT_log, "a", encoding="utf-8") as log_file:
        timestamp = datetime.now().isoformat(timespec="seconds")
        log_file.write(f"{timestamp} - {run_status}\n")

    if run_status != "success":
        print(f"CDFT 计算失败，请查看 '{CDFT_log}'")
        return "run_failed", None

    return "success", CDFT_txt


def get_excited_data(sub_dir):
    if not Path(sub_dir).exists():
        msg = f"目录 '{sub_dir}' 不存在！"
        print(f"错误: {msg}")
        return msg, None
    excited_log_dir = extract_last_task(sub_dir)
    excited_dir = Path(sub_dir / "excited.fchk")
    missing = [
        str(path.name) for path in (excited_log_dir, excited_dir) if not path.exists()
    ]
    if missing:
        msg = f"缺少 {', '.join(missing)}"
        print(f"错误: {sub_dir} {msg}")
        return msg, None
    input_text = f"""\
        {excited_dir}
        18
        1
        {excited_log_dir}
        1
        1
        2
        q
        """
    input_file = Path(sub_dir / "part2_input.txt")
    with open(input_file, "w", encoding="utf-8") as f:
        f.write(input_text)
    excited_status_log = Path(sub_dir / "excited_run.log")
    run_status = "failed"

    try:
        with open(input_file, "r", encoding="utf-8") as stdin_file, open(
            excited_status_log, "w", encoding="utf-8"
        ) as excited_stream:
            subprocess.run(
                ["Multiwfn_noGUI"],
                stdin=stdin_file,
                stdout=excited_stream,
                stderr=subprocess.STDOUT,
                cwd=sub_dir,
                text=True,
            )
    except Exception as e:
        msg = f"运行异常: {e}"
        print(f"错误: {msg}")
        return msg, None

    with open(excited_status_log, "r", encoding="utf-8") as f:
        content = f.read()
        if "Calculation took up wall clock time" in content:
            run_status = "success"

    with open(excited_status_log, "a", encoding="utf-8") as log_file:
        timestamp = datetime.now().isoformat(timespec="seconds")
        log_file.write(f"{timestamp} - {run_status}\n")
    if run_status != "success":
        print(f"激发态计算失败，请查看 '{excited_status_log}'")
        return "run_failed", None

    print(f"激发态输出位于 '{excited_status_log}'")
    return "success", excited_status_log


def data_to_csv(sub_dir, work_dir, CDFT_txt, excited_status_log):
    output_csv = Path(work_dir / "results.csv")
    excited_log_dir = Path(sub_dir / "excited.log")
    results = {}
    results["fchknam"] = str(sub_dir)

    if CDFT_txt is None or excited_status_log is None:
        print(f"跳过 {sub_dir}，输入文件缺失，未写入 CSV。")
        return

    if not Path(CDFT_txt).exists():
        print(f"跳过 {sub_dir}，未找到 {CDFT_txt}。")
        return

    if not Path(excited_status_log).exists():
        print(f"跳过 {sub_dir}，未找到 {excited_status_log}。")
        return

    if not excited_log_dir.exists():
        print(f"跳过 {sub_dir}，未找到 {excited_log_dir}。")
        return

    CDFT_txt_data = {
        "E_neutral_Hartree": r"E\(N\):\s+([-+]?\d*\.\d+)\s+Hartree",
        "E_anion_Hartree": r"E\(N\+1\):\s+([-+]?\d*\.\d+)\s+Hartree",
        "E_cation_Hartree": r"E\(N-1\):\s+([-+]?\d*\.\d+)\s+Hartree",
        "mulliken negativity": r"Mulliken electronegativity:\s+([-+]?\d*\.\d+)\s+Hartree",
        # "mulliken negativity_eV": r"Mulliken electronegativity:.*?([-+]?\d*\.\d+)\s+eV",
        "elec_index": r"Electrophilicity index:\s+([-+]?\d*\.\d+)\s+Hartree",
        # "elec_index_eV": r"Electrophilicity index:.*?([-+]?\d*\.\d+)\s+eV",
        "hardness": r"Hardness \(=fundamental gap\):\s+([-+]?\d*\.\d+)\s+Hartree",
        # "hardness_eV": r"Hardness \(=fundamental gap\):.*?([-+]?\d*\.\d+)\s+eV",
    }
    excited_txt_data = {
        "sr_index": r"Sr index \(integral of Sr function\):\s+([-+]?\d*\.\d+)\s",
        "d_index": r"D index:\s+([-+]?\d*\.\d+)",
        "t_index": r"t index:\s+([-+]?\d*\.\d+)",
        "hdi": r"Integral of hole:\s+([-+]?\d*\.\d+)",
        "edi": r"Integral of electron:\s+([-+]?\d*\.\d+)",    
    }
    excited_log_data = {
        "ex_en":r"Excited State\s+1:\s+Singlet-A\s+([-+]?\d+\.\d+)\s+eV",
        "fosc":r"Excited State\s+1:\s+Singlet-A.*?f=([-+]?\d+\.\d+)",
    }
    def to_float(value):
        try:
            return float(value)
        except (ValueError, TypeError):
            return None

    try:
        with open(excited_log_dir, "r", encoding="utf-8") as f:
            excited_log_content = f.read()

        for key, pattern in excited_log_data.items():
            match = re.search(pattern, excited_log_content)
            if match:
                results[key] = to_float(match.group(1))
            else:
                results[key] = None
    except Exception as e:
        print(f"Error reading excited.log file: {e}")
        for key in excited_log_data.keys():
            results[key] = None
            
    try:
        with open(CDFT_txt, "r", encoding="utf-8") as f:
            cdft_content = f.read()

        for key, pattern in CDFT_txt_data.items():
            match = re.search(pattern, cdft_content)
            if match:
                results[key] = to_float(match.group(1))
            else:
                results[key] = None
    except Exception as e:
        print(f"Error reading CDFT file: {e}")
        for key in CDFT_txt_data.keys():
            results[key] = None

    try:
        with open(excited_status_log, "r", encoding="utf-8") as f:
            excited_content = f.read()

        for key, pattern in excited_txt_data.items():
            match = re.search(pattern, excited_content)
            if match:
                results[key] = to_float(match.group(1))
            else:
                results[key] = None
    except Exception as e:
        print(f"Error reading excited log file: {e}")
        for key in excited_txt_data.keys():
            results[key] = None
    # 计算垂直电离能 (垂直电离能 = E(N-1) - E(N))
    if (
        results.get("E_cation_Hartree") is not None
        and results.get("E_neutral_Hartree") is not None
    ):
        results["vip"] = (
            results["E_cation_Hartree"] - results["E_neutral_Hartree"]
        )
    else:
        results["vip"] = None

    # 计算垂直电子亲和能 (垂直电子亲和能 = E(N) - E(N+1))
    if (
        results.get("E_neutral_Hartree") is not None
        and results.get("E_anion_Hartree") is not None
    ):
        results["vea"] = (
            results["E_neutral_Hartree"] - results["E_anion_Hartree"]
        )
    else:
        results["vea"] = None

    # 计算平均色散 (H index) = (λ_h + λ_e) / 2
    if (
        results.get("hdi") is not None
        and results.get("edi") is not None
    ):
        results["h_index"] = (
            results["hdi"] + results["edi"]
        ) / 2
    else:
        results["h_index"] = None

    # 计算色散差值 = λ_e - λ_h
    if (
        results.get("edi") is not None
        and results.get("hdi") is not None
    ):
        results["delta_sigma"] = (
            results["edi"] - results["hdi"]
        )
    else:
        results["delta_sigma"] = None

    ordered_data = {
        "fchknam": results["fchknam"],
        "ex_energy": results.get("ex_en"),
        "fosc": results.get("fosc"),
        "sr_index": results.get("sr_index"),
        "d_index": results.get("d_index"),
        "delta_sigma": results.get("delta_sigma"),
        "h_index": results.get("h_index"),
        "t_index": results.get("t_index"),
        "hdi": results.get("hdi"),
        "edi": results.get("edi"),
        "vip": results.get("vip"),
        "vea": results.get("vea"),
        "mulliken negativity": results.get("mulliken negativity"),
        "hardness": results.get("hardness"),
        "elec_index": results.get("elec_index"),
    }

    # 6. 写入CSV文件
    try:
        file_exists = output_csv.exists()
        
        # 定义CSV列的顺序
        field_order = [
            "fchknam", "ex_energy", "fosc", "sr_index", "d_index", 
            "delta_sigma", "h_index", "t_index", "hdi", "edi", 
            "vip", "vea", "mulliken negativity", "hardness", "elec_index"
        ]
        
        with open(output_csv, "a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=field_order)
            
            if not file_exists:
                writer.writeheader()
            
            # 确保所有字段都存在，用"N/A"填充缺失值
            row = {}
            for field in field_order:
                value = ordered_data.get(field)
                if value is None:
                    row[field] = "N/A"
                else:
                    row[field] = value
            
            writer.writerow(row)
        
        print(f"Results saved to {output_csv}")
        
    except Exception as e:
        print(f"Error writing CSV file: {e}")

    return


def load_existing_entries(csv_path):
    existing = set()
    if not Path(csv_path).exists():
        return existing

    try:
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = row.get("fchknam")
                if key:
                    existing.add(key)
    except Exception as e:
        print(f"读取 {csv_path} 失败: {e}")
    return existing


def write_pending_jobs(work_dir, pending_jobs):
    pending_file = Path(work_dir / "pending_jobs.csv")
    fieldnames = [
        "pair",
        "fragment",
        "sub_dir",
        "cdft_status",
        "excited_status",
        "updated_at",
    ]
    try:
        with open(pending_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for job in pending_jobs:
                writer.writerow(job)
        print(f"待补算列表写入 {pending_file}")
    except Exception as e:
        print(f"写入 {pending_file} 失败: {e}")


def main():
    work_dir = Path(config.work_dir)
    pair_dirs = work_dir.glob("pair*")
    results_csv = Path(work_dir / "results.csv")
    completed_entries = load_existing_entries(results_csv)
    pending_jobs = []

    for pair in pair_dirs:
        for sub_dir_name in ["amine", "aldehyde"]:
            sub_dir = pair / sub_dir_name
            cdft_status, CDFT_txt = get_CDFT_data(sub_dir)
            excited_status, excited_status_log = get_excited_data(sub_dir)

            if cdft_status == "success" and excited_status == "success":
                key = str(sub_dir)
                if key in completed_entries:
                    print(f"{sub_dir} 已存在于 results.csv，跳过写入。")
                    continue
                data_to_csv(sub_dir, work_dir, CDFT_txt, excited_status_log)
                completed_entries.add(key)
            else:
                pending_jobs.append(
                    {
                        "pair": pair.name,
                        "fragment": sub_dir_name,
                        "sub_dir": str(sub_dir),
                        "cdft_status": cdft_status,
                        "excited_status": excited_status,
                        "updated_at": datetime.now().isoformat(timespec="seconds"),
                    }
                )

    write_pending_jobs(work_dir, pending_jobs)
    return


if __name__ == "__main__":
    main()
