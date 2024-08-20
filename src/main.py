from databases import Database
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, inspect, select
from fastapi import FastAPI, HTTPException
from rdkit import Chem
import os
from dotenv import load_dotenv

app = FastAPI()

load_dotenv(".env")

DB_URL = os.getenv("DB_URL")

if DB_URL is None:
    raise ValueError("DATABASE_URL is not set in the environment variables")

# Debug print to check if the DB_URL is correct
print(f"DB_URL: {DB_URL}")

database = Database(DB_URL)
metadata = MetaData()

molecules = Table(
    "molecules",
    metadata,
    Column("identifier", Integer, primary_key=True),
    Column("smiles", String, nullable=False)
)

engine = create_engine(DB_URL)
inspector = inspect(engine)
if not inspector.has_table('molecules'):
    metadata.create_all(engine)

async def substructure_search(mols, mol):
    substructure = Chem.MolFromSmiles(mol)
    if substructure is None:
        raise ValueError(f'Invalid substructure SMILES: {mol}')
    substructure_num_atoms = substructure.GetNumAtoms()
    matches = []
    for molecule in mols:
        object_mol = Chem.MolFromSmiles(molecule)
        if object_mol is None:
            raise ValueError("Invalid Molecule!")
        if object_mol.GetNumAtoms() < substructure_num_atoms:
            continue
        if object_mol.HasSubstructMatch(substructure):
            matches.append(molecule)
    return matches

@app.on_event("startup")
async def startup():
    await database.connect()

@app.on_event("shutdown")
async def shutdown():
    await database.disconnect()

@app.get("/")
async def get_server():
    return {"server_id": os.getenv("SERVER_ID", "1")}

@app.post("/add")
async def add_molecule(id: int, smiles: str):
    query = molecules.insert().values(identifier=id, smiles=smiles)
    await database.execute(query)
    return {"message": f"Molecule '{id}' added successfully."}

@app.get("/get")
async def get_molecule(id: int):
    query = select(molecules.c.smiles).where(molecules.c.identifier == id)
    result = await database.fetch_one(query)
    if result:
        return result["smiles"]
    else:
        raise HTTPException(status_code=404, detail="Molecule not found")

@app.put("/update")
async def update_molecule(id: int, smiles: str):
    query = molecules.update().where(molecules.c.identifier == id).values(smiles=smiles)
    await database.execute(query)
    return {"message": f"Molecule '{id}' updated successfully."}

@app.delete("/del")
async def delete_molecule(id: int):
    query = molecules.delete().where(molecules.c.identifier == id)
    await database.execute(query)
    return {"message": f"Molecule '{id}' deleted successfully."}

@app.get("/getall")
async def list_all_molecules():
    query = select(molecules)
    rows = await database.fetch_all(query)
    return rows

@app.get("/subsearch")
async def sub_search(substructure: str):
    query = select(molecules.c.smiles)
    rows = await database.fetch_all(query)
    molecules_list = [row["smiles"] for row in rows]
    matches = await substructure_search(molecules_list, substructure)
    return matches