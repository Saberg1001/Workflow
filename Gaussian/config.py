from pathlib import Path

work_dir = Path(__file__).parent
script_file = work_dir /"tmp"/ "submit.sh"
tmp_file = work_dir /"tmp"/ "template.gjf"
input_log_file = work_dir /"log"/ "failed_jobs.log"
smiles_data_file = work_dir / "smiles_data.py"
generate_dir = work_dir / "data"/"input"
ana_log_dir = work_dir / "log"/ "analysis.log"
part1_dir = work_dir / "part1.txt"
part2_dir = work_dir / "part2.txt"
       