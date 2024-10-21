from sqlalchemy.orm import Session
from sqlalchemy.exc import SQLAlchemyError
from fastapi import FastAPI, HTTPException, Depends
from celery.result import AsyncResult
from celery_worker import celery
from tasks import substructure_search_task
from utils import logger
from schemas import MoleculeCreate, MoleculeUpdate, MoleculeSchema 
from models import Molecule as MoleculeModel
from crud import MoleculeDAO
from typing import List
from database import get_db
from redis_utils import get_cached_result
from rdkit import Chem
import os


app = FastAPI()


@app.get("/")
async def get_server():
    server_id = os.getenv("SERVER_ID", "1")
    logger.info(f"Server ID: {server_id}")
    return {"server_id": server_id}


def is_valid_smiles(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

@app.post("/add", response_model=MoleculeSchema)
async def add_molecule(molecule: MoleculeCreate, db: Session = Depends(get_db)):
    if not is_valid_smiles(molecule.smiles):
        logger.warning(f"Invalid SMILES string: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    molecule_dao = MoleculeDAO(db)
    db_molecule = molecule_dao.create_molecule(molecule)
    logger.info(f"Molecule '{db_molecule.id}' added successfully.")
    return db_molecule
    # return {"message": f"Molecule '{molecule.id}' added successfully."}


@app.get("/get", response_model=MoleculeSchema)
async def get_molecule(id: int, db: Session = Depends(get_db)):
    molecule_dao = MoleculeDAO(db)
    db_molecule = molecule_dao.get_molecule(id)
    if db_molecule is None:
        logger.warning(f"Molecule '{id}' not found.")
        raise HTTPException(status_code=404, detail="Molecule not found")
    logger.info(f"Molecule '{id}' found.")
    return db_molecule


@app.put("/update", response_model=MoleculeSchema)
async def update_molecule(id: int, molecule: MoleculeUpdate, db: Session = Depends(get_db)):
    molecule_dao = MoleculeDAO(db)
    db_molecule = molecule_dao.get_molecule(id)
    if db_molecule is None:
        logger.warning(f"Molecule '{id}' not found.")
        raise HTTPException(status_code=404, detail="Molecule not found")

    if molecule.smiles and not is_valid_smiles(molecule.smiles):
        logger.warning(f"Invalid SMILES string: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
       
    db_molecule = molecule_dao.update_molecule(db_molecule, molecule)
    logger.info(f"Molecule '{id}' updated successfully.")
    return db_molecule
    # return {"message": f"Molecule '{molecule.id}' updated successfully."}


@app.delete("/del")
async def delete_molecule(id: int, db: Session = Depends(get_db)):
    molecule_dao = MoleculeDAO(db)
    db_molecule = molecule_dao.get_molecule(id)
    if db_molecule is None:
        logger.warning(f"Molecule '{id}' not found.")
        raise HTTPException(status_code=404, detail="Molecule not found")
    molecule_dao.delete_molecule(db_molecule)
    logger.info(f"Molecule '{id}' deleted successfully.")
    return {"message": f"Molecule '{id}' deleted successfully."}


@app.get("/getall", response_model=List[MoleculeSchema])
def list_all_molecules(limit: int = 100, offset: int = 0, db: Session = Depends(get_db)):
    try:
        molecules = db.query(MoleculeModel).offset(offset).limit(limit).all()
        logger.info(f"Fetched {len(molecules)} molecules from the database.")
        return molecules
    except SQLAlchemyError as e:
        logger.error(f"An error occurred: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")


@app.get("/subsearch")
async def sub_search(substructure: str, db: Session = Depends(get_db)):
    cache_key = f"subsearch:{substructure}"
    try:
        cached_result = get_cached_result(cache_key)

        if cached_result is not None:
            logger.info(f"Returning cached result for substructure: {substructure}")
            return {"source": "cache", "data": cached_result}

        molecules = db.query(MoleculeModel.smiles).all()
        molecules_list = [m.smiles for m in molecules]
        task = substructure_search_task.apply_async(args=[substructure, molecules_list])    

        logger.info(f"Substructure search task initiated with task ID: {task.id}")
        return {"task_id": task.id, "status": "PENDING"}

    except Exception as e:
        logger.error(f"An error occurred during substructure search: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")


@app.get("/subsearch/{task_id}")
async def get_task_result(task_id: str):
    task_result = AsyncResult(task_id, app=celery)
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {"task_id": task_id, "status": "Task completed", "result": task_result.result}
    else:
        return {"task_id": task_id, "status": task_result.state}
