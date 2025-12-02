from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from smiles_data import pairs


def smiles_to_gjf(smile, name, work_dir):
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
    gjf_dir = Path(work_dir / f"{name}.gjf")
    with open(gjf_dir, "w") as f:
        f.write("#p opt b3lyp/6-31g(d) \n\n")
        f.write("title\n")
        f.write("\n0 1\n")
        for atom in mol_h.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            f.write(f"{symbol} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")
    return gjf_dir


def repalce_coordinates(work_dir, gjf_dir, smile, name, temp_file, submit_script):
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
    with open(temp_file, "r",encoding='utf-8') as f:
        content = f.read()
        new_content = content.replace("atomic coordinates", atom_coord)
    input_name = f"{name}-input.gjf"
    input_dir = work_dir / input_name
    with open(input_dir, "w") as f:
        f.write(new_content)
    
    with open(submit_script, "r",encoding='utf-8') as f:
        script_content = f.read()
        new_script_content = script_content.replace("jobname", f"{gjf_dir}")
    submit_dir = work_dir / f"{name}.sh"
    with open(submit_dir, "w") as f:
        f.write(new_script_content)
    return 



def gjf_to_link_input(tmp_file="template.gjf", submit_script="submit.sh"):
    for pair in pairs:
        if pair["name"] == "pair01":
            name = pair["name"]
            work_dir = Path(__file__).parent/name
            aldehyde = pair["aldehyde"]
            amine = pair["amine"]
            if not Path(work_dir).exists():
                Path(work_dir).mkdir(parents=True, exist_ok=True)
            aldehyde_work_dir = work_dir / "aldehyde"
            if not Path(aldehyde_work_dir).exists():
                Path(aldehyde_work_dir).mkdir(parents=True, exist_ok=True)
            aldehyde_gjf_dir = smiles_to_gjf(aldehyde,"aldehyde",aldehyde_work_dir)
            repalce_coordinates(aldehyde_work_dir,aldehyde_gjf_dir,aldehyde,"aldehyde", tmp_file,submit_script)
            
            amine_work_dir = work_dir / "amine"
            if not Path(amine_work_dir).exists():
                Path(amine_work_dir).mkdir(parents=True, exist_ok=True)
                 
            amine_gjf_dir = smiles_to_gjf(amine,"amine", amine_work_dir)    
            repalce_coordinates(amine_work_dir,amine_gjf_dir,amine,"aldehyde", tmp_file,submit_script)
        return


if __name__ == "__main__":
    gjf_to_link_input()
