from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from datetime import datetime
import config
from smiles_data import pairs


def smiles_to_gjf(smile, name, work_dir):
    gjf_dir = Path(work_dir / f"{name}.gjf")
    if gjf_dir.exists():
        print("GJF file already exists, skipping generation.")
    else:
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        mol_h = Chem.AddHs(mol)

        state = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
        if state != 0:
            raise ValueError("3D embedding failed")
        try:
            AllChem.UFFOptimizeMolecule(mol_h)
        except:
            print("UFF optimization failed, proceeding without optimization.")

        conf = mol_h.GetConformer()
        
        with open(gjf_dir, "w") as f:
            f.write("#p opt b3lyp/6-31g(d) \n\n")
            f.write("title\n")
            f.write("\n0 1\n")
            for atom in mol_h.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                symbol = atom.GetSymbol()
                f.write(f"{symbol} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")
    return gjf_dir


def log_failure(log_path, pair_name, species, smile, error_message):
    timestamp = datetime.now().isoformat(timespec="seconds")
    with open(log_path, "a", encoding="utf-8") as log:
        log.write(f"{timestamp}\t{pair_name}\t{species}\t{smile}\t{error_message}\n")


def repalce_coordinates(work_dir, gjf_dir, smile, name, tmp_file, submit_script):
    coord = []
    with open(gjf_dir, "r",encoding='utf-8') as f:
        lines = f.readlines()
        for line in lines:
            if line.strip() and len(line.split()) == 4:
                try:
                    float(line.split()[1])
                    coord.append(line)
                except:
                    pass
    atom_coord = "".join(coord)
    with open(tmp_file, "r",encoding='utf-8') as f:
        content = f.read()
        new_content = content.replace("atomic coordinates", atom_coord)
    input_name = f"{name}-input.gjf"
    input_dir = work_dir / input_name
    with open(input_dir, "w") as f:
        f.write(new_content)
    
    with open(submit_script, "r",encoding='utf-8') as f:
        script_content = f.read()
        new_script_content = script_content.replace("jobname", input_name)
    submit_dir = work_dir / "submit.sh"
    with open(submit_dir, "w") as f:
        f.write(new_script_content)
    return 



def gjf_to_link_input():
    work_dir = config.work_dir
    if not Path(work_dir).exists():
        Path(work_dir).mkdir(parents=True, exist_ok=True)
    tmp_file = config.tmp_file
    submit_script = config.script_file
    log_file = config.input_log_file
    generate_dir = config.generate_dir
    failures = []
    print(tmp_file,submit_script,log_file,work_dir)
    for pair in pairs:
        #if pair["name"] == "pair27":
            name = pair["name"]
            pair_dir = generate_dir/name
            aldehyde = pair["aldehyde"]
            amine = pair["amine"]
            if not Path(pair_dir).exists():
                Path(pair_dir).mkdir(parents=True, exist_ok=True)

            species_configs = [
                ("aldehyde", aldehyde, pair_dir / "aldehyde"),
                ("amine", amine, pair_dir / "amine"),
            ]

            for species, smile, target_dir in species_configs:
                if not target_dir.exists():
                    target_dir.mkdir(parents=True, exist_ok=True)
                try:
                    gjf_dir = smiles_to_gjf(smile, species, target_dir)
                    repalce_coordinates(target_dir, gjf_dir, smile, species, tmp_file, submit_script)
                except Exception as exc:
                    message = str(exc)
                    print(f"{name} 的 {species} 生成失败，已跳过：{message}")
                    log_failure(log_file, name, species, smile, message)
                    failures.append((name, species, smile, message))
    if failures:
        print(f"共 {len(failures)} 个分子需要手动生成，详情见 {log_file}")
    return

if __name__ == "__main__":
    gjf_to_link_input()
