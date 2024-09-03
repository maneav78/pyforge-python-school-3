import logging
from databases import Database
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, inspect, select
from fastapi import FastAPI, HTTPException, Depends
from contextlib import asynccontextmanager
from rdkit import Chem
import os
from dotenv import load_dotenv
from pydantic import BaseModel

app = FastAPI()

load_dotenv(".env")

DB_URL = os.getenv("DB_URL")

if DB_URL is None:
    raise ValueError("DATABASE_URL is not set in the environment variables")

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


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Molecule(BaseModel):
    id: int
    smiles: str


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


@asynccontextmanager
async def lifespan(app: FastAPI):
    await database.connect()
    logger.info("Database connected")
    yield
    await database.disconnect()
    logger.info("Database disconnected")

app.router.lifespan_context = lifespan


async def get_db() -> Database:
    return database


@app.get("/")
async def get_server():
    server_id = os.getenv("SERVER_ID", "1")
    logger.info(f"Server ID: {server_id}")
    return {"server_id": server_id}


@app.post("/add")
async def add_molecule(molecule: Molecule, db: Database = Depends(get_db)):
    query = molecules.insert().values(identifier=molecule.id, smiles=molecule.smiles)
    await db.execute(query)
    logger.info(f"Molecule '{molecule.id}' added successfully.")
    return {"message": f"Molecule '{molecule.id}' added successfully."}


@app.get("/get")
async def get_molecule(id: int, db: Database = Depends(get_db)):
    query = select(molecules.c.smiles).where(molecules.c.identifier == id)
    result = await db.fetch_one(query)
    if result:
        logger.info(f"Molecule '{id}' found.")
        return result["smiles"]
    else:
        logger.warning(f"Molecule '{id}' not found.")
        raise HTTPException(status_code=404, detail="Molecule not found")


@app.put("/update")
async def update_molecule(molecule: Molecule, db: Database = Depends(get_db)):
    query = molecules.update().where(molecules.c.identifier == molecule.id).values(smiles=molecule.smiles)
    await db.execute(query)
    logger.info(f"Molecule '{molecule.id}' updated successfully.")
    return {"message": f"Molecule '{molecule.id}' updated successfully."}


@app.delete("/del")
async def delete_molecule(id: int, db: Database = Depends(get_db)):
    query = molecules.delete().where(molecules.c.identifier == id)
    await db.execute(query)
    logger.info(f"Molecule '{id}' deleted successfully.")
    return {"message": f"Molecule '{id}' deleted successfully."}


@app.get("/getall")
async def list_all_molecules(limit: int = 100, db: Database = Depends(get_db)):
    try:
        query = select(molecules)
        rows = await db.fetch_all(query)
        logger.info(f"Fetched {len(rows)} molecules from the database.")

        count = 0
        results = []
        for row in rows:
            if count >= limit:
                break
            results.append(row)
            count += 1
            logger.info(f"Yielding molecule {count} with ID: {row['identifier']}")

        return results
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")


@app.get("/subsearch")
async def sub_search(substructure: str, db: Database = Depends(get_db)):
    query = select(molecules.c.smiles)
    rows = await db.fetch_all(query)
    molecules_list = [row["smiles"] for row in rows]
    matches = await substructure_search(molecules_list, substructure)
    logger.info(f"Substructure search completed with {len(matches)} matches.")
    return matches
