from unittest.mock import AsyncMock, MagicMock, patch
from fastapi.testclient import TestClient
from src.main import app, get_db
from dotenv import load_dotenv
import os
import psycopg2
from psycopg2 import sql
import pytest
from sqlalchemy import create_engine, text
from sqlalchemy.exc import OperationalError
from databases import Database

load_dotenv(".env")

MAIN_DATABASE_URL = os.getenv("MAIN_DATABASE_URL")
TEST_DB_NAME = os.getenv("TEST_DB_NAME")
TEST_DATABASE_URL = os.getenv("TEST_DATABASE_URL")

with psycopg2.connect(MAIN_DATABASE_URL) as conn:
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL("SELECT 1 FROM pg_database WHERE datname = %s"),
            [TEST_DB_NAME]
        )
        exists = cur.fetchone()
        if not exists:
            cur.execute(f"CREATE DATABASE {TEST_DB_NAME}")

test_engine = create_engine(TEST_DATABASE_URL)
test_database = Database(TEST_DATABASE_URL)

@pytest.fixture
def db_connection():
    try:
        connection = test_engine.connect()
        yield connection
    finally:
        connection.close()

def test_postgres_connection(db_connection):
    try:
        result = db_connection.execute(text("SELECT 1"))
        assert result.fetchone() == (1,)
    except OperationalError as e:
        pytest.fail(f"Database connection failed: {e}")

client = TestClient(app)

def override_get_db(initial_fetch_one_value=None):
    db = MagicMock()
    db.fetch_one = AsyncMock(return_value=initial_fetch_one_value)
    db.execute = AsyncMock(return_value=None)
    return db

def setup_db_override(initial_fetch_one_value=None):
    db = override_get_db(initial_fetch_one_value)
    app.dependency_overrides[get_db] = lambda: db
    return db

def add_molecule(client, molecule_id, smiles):
    response = client.post("/add", json={"id": molecule_id, "smiles": smiles})
    assert response.status_code == 200
    assert response.json() == {"message": f"Molecule '{molecule_id}' added successfully."}

def test_add_molecule():
    setup_db_override()
    add_molecule(client, 1, "CCO")

def test_get_molecule():
    setup_db_override({"smiles": "CCO"})
    add_molecule(client, 1, "CCO")
    response = client.get("/get", params={"id": 1})
    assert response.status_code == 200
    assert response.json() == "CCO"

def test_update_molecule():
    db = setup_db_override({"smiles": "CCO"})
    add_molecule(client, 3, "CCO")
    response = client.put("/update", json={"id": 3, "smiles": "c1ccccc1"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule '3' updated successfully."}
    db.fetch_one.return_value = {"smiles": "c1ccccc1"}
    response = client.get("/get", params={"id": 3})
    assert response.status_code == 200
    assert response.json() == "c1ccccc1"

def test_delete_molecule():
    db = setup_db_override({"smiles": "CCO"})
    add_molecule(client, 4, "CCO")
    response = client.delete("/del", params={"id": 4})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule '4' deleted successfully."}
    db.fetch_one.return_value = None
    response = client.get("/get", params={"id": 4})
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}

def mock_fetch_all_return_values(db, values):
    async def mock_fetch_all(query):
        return values
    db.fetch_all = mock_fetch_all

def test_list_all_molecules():
    db = setup_db_override()
    molecules = [
        {"id": 1, "smiles": "CCO"},
        {"id": 2, "smiles": "C1=CC=CC=C1"},
        {"id": 3, "smiles": "NCCO"}
    ]
    add_molecule(client, 1, "CCO")
    add_molecule(client, 2, "C1=CC=CC=C1")
    add_molecule(client, 3, "NCCO")
    mock_fetch_all_return_values(db, molecules)
    response = client.get("/getall")
    assert response.status_code == 200
    assert response.json() == molecules

def test_list_all_molecules_empty():
    db = setup_db_override()
    mock_fetch_all_return_values(db, [])
    response = client.get("/getall")
    assert response.status_code == 200
    assert response.json() == []

@patch("src.main.substructure_search")
def test_sub_search(mock_substructure_search):
    db = setup_db_override()
    molecules = [
        {"smiles": "CCO"},
        {"smiles": "C1=CC=CC=C1"},
        {"smiles": "NCCO"}
    ]
    mock_fetch_all_return_values(db, molecules)
    mock_substructure_search.return_value = ["CCO", "NCCO"]
    response = client.get("/subsearch", params={"substructure": "CCO"})
    assert response.status_code == 200
    assert response.json() == ["CCO", "NCCO"]

def test_sub_search_no_matches():
    db = setup_db_override()
    molecules = [
        {"smiles": "CCO"},
        {"smiles": "C1=CC=CC=C1"},
        {"smiles": "NCCO"}
    ]
    mock_fetch_all_return_values(db, molecules)
    with patch("src.main.substructure_search", return_value=[]):
        response = client.get("/subsearch", params={"substructure": "XYZ"})
        assert response.status_code == 200
        assert response.json() == []
