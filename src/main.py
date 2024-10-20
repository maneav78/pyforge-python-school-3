from databases import Database
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, inspect, select
from fastapi import FastAPI, HTTPException, Depends
from contextlib import asynccontextmanager
import os
from dotenv import load_dotenv
from pydantic import BaseModel
import redis
import json
from celery.result import AsyncResult
from celery_worker import celery
from tasks import substructure_search_task
from utils import logger


app = FastAPI()

redis_client = redis.Redis(host='redis', port=6379, db=0)

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


class Molecule(BaseModel):
    id: int
    smiles: str


# redis functions
def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: dict, expiration: int = 60):
    redis_client.setex(key, expiration, json.dumps(value))


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
    cache_key = f"subsearch:{substructure}"
    try:
        cached_result = get_cached_result(cache_key)

        if cached_result is not None:
            logger.info(f"Returning cached result for substructure: {substructure}")
            return {"source": "cache", "data": cached_result}

        query = select(molecules.c.smiles)
        rows = await db.fetch_all(query)
        molecules_list = [row["smiles"] for row in rows]
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
