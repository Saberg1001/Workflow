import glob
from pathlib import Path
import subprocess

pair_dirs = Path("./").glob("pair*")

for pair in pair_dirs:
    for sub_dir_name in ["amine", "aldehyde"]:
        sub_dir = pair / sub_dir_name
        submit_script = sub_dir / "submit.sh"
        if not submit_script.is_file():
            print(f"Submitting in {sub_dir.parent}")
            
            subprocess.run(
                ["sbatch", "submit.sh"],
                cwd=str(sub_dir.parent),
                check=True
            )