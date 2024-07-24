from rdkit import Chem
from rdkit.Chem import Draw

def substructure_search(mols, mol):
    substructure = Chem.MolFromSmiles(mol)

    if substructure is None:
        raise ValueError(f'Invalid substructure SMILES: {mol}')
    
    substructure_num_atoms = substructure.GetNumAtoms()

    matches = []
    for molecule in mols:
        object_mol = Chem.MolFromSmiles(molecule)
        if object_mol is None:
            raise ValueError(f'Invalid Molecule!')

        if object_mol.GetNumAtoms() < substructure_num_atoms:
            continue 
        
        if object_mol.HasSubstructMatch(substructure):
           matches.append(molecule)

    return matches

# testing
if name == "main":
    assert substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1") == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
    
    assert substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "invalid_smiles") == []
    
    print("All tests passed!")