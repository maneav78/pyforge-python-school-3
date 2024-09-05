from unittest.mock import AsyncMock, MagicMock, patch
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, text, inspect, MetaData
from sqlalchemy.exc import OperationalError
from databases import Database
from dotenv import load_dotenv
import os
import psycopg2
from psycopg2 import sql
import pytest
from httpx import AsyncClient
from src.main import app, get_db

load_dotenv(".env")

MAIN_DATABASE_URL = os.getenv("MAIN_DATABASE_URL")
TEST_DB_NAME = os.getenv("TEST_DB_NAME")
TEST_DATABASE_URL = os.getenv("TEST_DATABASE_URL")

exists = None


def check_database_exists(conn, db_name):
    global exists
    cur = conn.cursor()
    try:
        cur.execute(
            sql.SQL("SELECT 1 FROM pg_database WHERE datname = %s"),
            [db_name]
        )
        exists = cur.fetchone()
    finally:
        cur.close()


def create_database(db_name):
    conn = psycopg2.connect(MAIN_DATABASE_URL)
    try:
        conn.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)
        cur = conn.cursor()
        try:
            cur.execute(sql.SQL("CREATE DATABASE {}").format(sql.Identifier(db_name)))
        finally:
            cur.close()
    finally:
        conn.close()


conn = psycopg2.connect(MAIN_DATABASE_URL)
try:
    conn.autocommit = True
    check_database_exists(conn, TEST_DB_NAME)
finally:
    conn.close()


test_engine = create_engine(TEST_DATABASE_URL)
test_database = Database(TEST_DATABASE_URL)
metadata = MetaData()

if not inspect(test_engine).has_table("molecules"):
    metadata.create_all(test_engine, checkfirst=True)


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


@pytest.fixture()
def client():
    return TestClient(app)  # Create a TestClient instance with your FastAPI app


def override_get_db(initial_fetch_one_value=None):
    db = MagicMock()
    db.fetch_one = AsyncMock(return_value=initial_fetch_one_value)
    db.execute = AsyncMock(return_value=None)
    db.fetch_all = AsyncMock(return_value=[])
    return db


def setup_db_override(initial_fetch_one_value=None):
    db = override_get_db(initial_fetch_one_value)
    app.dependency_overrides[get_db] = lambda: db
    return db


def add_molecule(client, molecule_id, smiles):
    response = client.post("/add", json={"id": molecule_id, "smiles": smiles})
    assert response.status_code == 200
    assert response.json() == {"message": f"Molecule '{molecule_id}' added successfully."}


def test_add_molecule(client):
    setup_db_override()
    add_molecule(client, 1, "CCO")


def test_get_molecule(client):
    setup_db_override({"smiles": "CCO"})
    add_molecule(client, 1, "CCO")
    response = client.get("/get", params={"id": 1})
    assert response.status_code == 200
    assert response.json() == "CCO"


def test_update_molecule(client):
    db = setup_db_override({"smiles": "CCO"})
    add_molecule(client, 3, "CCO")
    response = client.put("/update", json={"id": 3, "smiles": "c1ccccc1"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule '3' updated successfully."}
    db.fetch_one.return_value = {"smiles": "c1ccccc1"}
    response = client.get("/get", params={"id": 3})
    assert response.status_code == 200
    assert response.json() == "c1ccccc1"


def test_delete_molecule(client):
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


def test_list_all_molecules(client):
    db = setup_db_override()
    molecules = [
        {"identifier": 1, "smiles": "CCO"},
        {"identifier": 2, "smiles": "C1=CC=CC=C1"},
        {"identifier": 3, "smiles": "NCCO"}
    ]
    add_molecule(client, 1, "CCO")
    add_molecule(client, 2, "C1=CC=CC=C1")
    add_molecule(client, 3, "NCCO")
    mock_fetch_all_return_values(db, molecules)
    response = client.get("/getall")
    assert response.status_code == 200
    assert response.json() == molecules


def test_list_all_molecules_empty(client):
    db = setup_db_override()
    mock_fetch_all_return_values(db, [])
    response = client.get("/getall")
    assert response.status_code == 200
    assert response.json() == []


@pytest.fixture
def mock_cache():
    with patch("src.main.get_cached_result") as mock_get_cached_result, \
         patch("src.main.set_cache") as mock_set_cache:
        yield mock_get_cached_result, mock_set_cache


@pytest.mark.asyncio
@patch("src.main.substructure_search", new_callable=AsyncMock)
async def test_sub_search_cache_hit(mock_substructure_search):
    async with AsyncClient(app=app, base_url="http://test") as client:
        substructure = "CCO"
        mock_substructure_search.return_value = ["CCO"]

        with patch("src.main.get_cached_result", return_value=["CCO"]), \
             patch("src.main.set_cache") as mock_set_cache:
            response = await client.get("/subsearch", params={"substructure": substructure})

            assert response.status_code == 200
            assert response.json() == {"source": "cache", "data": ["CCO"]}
            mock_substructure_search.assert_not_called()


@pytest.mark.asyncio
async def test_sub_search_cache_miss(client, mock_cache):
    db = setup_db_override()
    mock_get_cached_result, mock_set_cache = mock_cache
    substructure = "CCO"
    mock_get_cached_result.return_value = None

    mock_rows = [{"smiles": "CCO"}, {"smiles": "C1=CC=CC=C1"}]
    db.fetch_all = AsyncMock(return_value=mock_rows)
    with patch("src.main.substructure_search", return_value=["CCO"]) as mock_substructure_search:
        response = client.get("/subsearch", params={"substructure": substructure})
        assert response.status_code == 200
        assert response.json() == {"source": "database", "data": ["CCO"]}
        mock_get_cached_result.assert_called_once_with(f"subsearch:{substructure}")
        mock_set_cache.assert_called_once_with(f"subsearch:{substructure}", ["CCO"], expiration=300)
        mock_substructure_search.assert_called_once_with(["CCO", "C1=CC=CC=C1"], substructure)
