import logging
from rdkit import Chem

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def substructure_search(mols, mol):
    logger.info("Starting substructure search")
    substructure = Chem.MolFromSmiles(mol)
    if substructure is None:
        logger.error(f'Invalid substructure SMILES: {mol}')
        raise ValueError(f'Invalid substructure SMILES: {mol}')
    substructure_num_atoms = substructure.GetNumAtoms()
    matches = []
    for molecule in mols:
        object_mol = Chem.MolFromSmiles(molecule)
        if object_mol is None:
            logger.error("Invalid Molecule!")
            raise ValueError("Invalid Molecule!")
        if object_mol.GetNumAtoms() < substructure_num_atoms:
            continue
        if object_mol.HasSubstructMatch(substructure):
            matches.append(molecule)
    logger.info(f"Substructure search complete, found {len(matches)} matches")
    return matches
