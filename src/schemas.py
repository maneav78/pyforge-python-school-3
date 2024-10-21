from pydantic import BaseModel


class MoleculeBase(BaseModel):
    smiles: str


class MoleculeCreate(MoleculeBase):
    pass


class MoleculeUpdate(MoleculeBase):
    pass


class MoleculeInDBBase(MoleculeBase):
    id: int

    class Config:
        from_attributes = True


class MoleculeSchema(MoleculeInDBBase):
    pass
