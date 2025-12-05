from pathlib import Path

work_dir = Path(__file__).parent
script_file = work_dir /"tmp"/ "submit.sh"
tmp_file = work_dir /"tmp"/ "template.gjf"
log_file = work_dir /"log"/ "failed_jobs.log"
smiles_data_file = work_dir / "smiles_data.py"
generate_dir = work_dir / "data"/"input"

print(f"work_dir: {work_dir}")
print(f"script_file: {script_file}")        