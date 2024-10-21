from sqlalchemy.orm import Session
import models
import schemas


class MoleculeDAO:
    def __init__(self, db: Session):
        self.db = db

    def create_molecule(self, molecule_in: schemas.MoleculeCreate):
        db_molecule = models.Molecule(smiles=molecule_in.smiles)
        self.db.add(db_molecule)
        self.db.commit()
        self.db.refresh(db_molecule)
        return db_molecule

    def get_molecule(self, molecule_id: int):
        return self.db.query(models.Molecule).filter(models.Molecule.id == molecule_id).first()

    def update_molecule(self, db_molecule: models.Molecule, molecule_in: schemas.MoleculeUpdate):
        db_molecule.smiles = molecule_in.smiles
        self.db.commit()
        self.db.refresh(db_molecule)
        return db_molecule

    def delete_molecule(self, db_molecule: models.Molecule):
        self.db.delete(db_molecule)
        self.db.commit()
