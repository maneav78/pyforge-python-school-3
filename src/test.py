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


def drop_existing_schema(engine):
    inspector = inspect(engine)
    # Drop all tables
    for table_name in inspector.get_table_names():
        engine.execute(f"DROP TABLE IF EXISTS {table_name} CASCADE")
    # Drop all sequences
    for seq_name in inspector.get_sequence_names():
        engine.execute(f"DROP SEQUENCE IF EXISTS {seq_name} CASCADE")


conn = psycopg2.connect(MAIN_DATABASE_URL)
try:
    conn.autocommit = True
    check_database_exists(conn, TEST_DB_NAME)
    if not exists:
        create_database(TEST_DB_NAME)
finally:
    conn.close()


test_engine = create_engine(TEST_DATABASE_URL)
test_database = Database(TEST_DATABASE_URL)
metadata = MetaData()


@pytest.fixture(scope="function", autouse=True)
def setup_and_teardown_db():
    drop_existing_schema(test_engine)
    metadata.create_all(test_engine, checkfirst=True)
    yield
    metadata.drop_all(test_engine)


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
    return TestClient(app)


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
async def test_sub_search_cache_hit(mock_cache):
    mock_get_cached_result, _ = mock_cache
    mock_get_cached_result.return_value = ["CCO"]
    async with AsyncClient(app=app, base_url="http://test") as client:
        substructure = "CCO"
        response = await client.get("/subsearch", params={"substructure": substructure})

        assert response.status_code == 200
        assert response.json() == {"source": "cache", "data": ["CCO"]}
        mock_get_cached_result.assert_called_once_with(f"subsearch:{substructure}")


@pytest.fixture
def mock_substructure_search_task():
    with patch("src.main.substructure_search_task") as mock:
        yield mock


@pytest.fixture
def mock_get_db():
    with patch("src.main.get_db") as mock:
        yield mock


def test_sub_search_without_cache(client, mock_cache, mock_substructure_search_task, mock_get_db):
    mock_get_cached_result, _ = mock_cache
    mock_get_cached_result.return_value = None

    mock_substructure_search_task.apply_async.return_value.id = "test_task_id"

    mock_db = mock_get_db.return_value
    mock_db.fetch_all.return_value = [
        {"identifier": 1, "smiles": "CCO"},
        {"identifier": 2, "smiles": "C1=CC=CC=C1"},
        {"identifier": 3, "smiles": "NCCO"}
    ]

    response = client.get("/subsearch", params={"substructure": "CCO"})

    assert response.status_code == 200
    assert response.json() == {"task_id": "test_task_id", "status": "PENDING"}
    assert mock_substructure_search_task.apply_async.called


@pytest.mark.asyncio
@patch('src.main.AsyncResult')
async def test_get_task_result_success(mock_async_result):
    task_id = 'task-id'

    mock_result = mock_async_result.return_value
    mock_result.state = 'SUCCESS'
    mock_result.result = {'data': 'some result'}

    async with AsyncClient(app=app, base_url="http://test") as client:
        response = await client.get(f"/subsearch/{task_id}")

    assert response.status_code == 200
    assert response.json() == {"task_id": task_id, "status": "Task completed", "result": {'data': 'some result'}}


@pytest.mark.asyncio
async def test_get_task_result_failed():
    task_id = 'task-id'

    mock_result = MagicMock()
    mock_result.state = 'FAILURE'
    mock_result.result = None

    with patch('src.main.AsyncResult', return_value=mock_result):
        async with AsyncClient(app=app, base_url="http://test") as client:
            response = await client.get(f"/subsearch/{task_id}")

    assert response.status_code == 200
    assert response.json() == {"task_id": task_id, "status": 'FAILURE'}


@pytest.mark.asyncio
@patch('src.main.AsyncResult')
async def test_get_task_result_pending(mock_async_result):
    task_id = 'task-id'

    mock_result = mock_async_result.return_value
    mock_result.state = 'PENDING'
    mock_result.result = None

    async with AsyncClient(app=app, base_url="http://test") as client:
        response = await client.get(f"/subsearch/{task_id}")

    assert response.status_code == 200
    assert response.json() == {"task_id": task_id, "status": "Task is still processing"}
