import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from main import app, get_db
from models import Base
from unittest.mock import patch, MagicMock
import os

SQLALCHEMY_TEST_DATABASE_URL = "sqlite:///./test.db"

engine = create_engine(
    SQLALCHEMY_TEST_DATABASE_URL, connect_args={"check_same_thread": False}
)
TestingSessionLocal = sessionmaker(bind=engine)

@pytest.fixture(scope="session", autouse=True)
def create_test_database():
    Base.metadata.create_all(bind=engine)
    yield
    Base.metadata.drop_all(bind=engine)

@pytest.fixture(scope="function")
def db_session():
    connection = engine.connect()
    transaction = connection.begin()
    session = TestingSessionLocal(bind=connection)

    # Override the get_db dependency to use the test session
    app.dependency_overrides[get_db] = lambda: session

    yield session

    session.close()
    transaction.rollback()
    connection.close()

client = TestClient(app)
def test_add_molecule():
    response = client.post("/add", json={"smiles": "CCO"})
    assert response.status_code == 200
    data = response.json()
    assert data["smiles"] == "CCO"
    assert "id" in data

def test_get_molecule():
    response = client.post("/add", json={"smiles": "C1=CC=CC=C1"})
    data = response.json()
    molecule_id = data["id"]

    response = client.get(f"/get?id={molecule_id}")
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == molecule_id
    assert data["smiles"] == "C1=CC=CC=C1"

def test_update_molecule():
    response = client.post("/add", json={"smiles": "CCO"})
    data = response.json()
    molecule_id = data["id"]

    response = client.put(f"/update?id={molecule_id}", json={"smiles": "C1=CC=CC=C1"})
    assert response.status_code == 200
    data = response.json()
    assert data["smiles"] == "C1=CC=CC=C1"

    response = client.get(f"/get?id={molecule_id}")
    data = response.json()
    assert data["smiles"] == "C1=CC=CC=C1"

def test_delete_molecule():
    response = client.post("/add", json={"smiles": "CCO"})
    data = response.json()
    molecule_id = data["id"]

    response = client.delete(f"/del?id={molecule_id}")
    assert response.status_code == 200
    data = response.json()
    assert data == {"message": f"Molecule '{molecule_id}' deleted successfully."}

    response = client.get(f"/get?id={molecule_id}")
    assert response.status_code == 404

def test_list_all_molecules():
    for smiles in ["CCO", "C1=CC=CC=C1", "NCCO"]:
        client.post("/add", json={"smiles": smiles})

    response = client.get("/getall")
    data = response.json()
    assert len(data) == 3  

    response = client.get("/getall?limit=2")
    data = response.json()
    assert len(data) == 2

    response = client.get("/getall?limit=1&offset=1")
    data = response.json()
    assert len(data) == 1
    assert data[0]["smiles"] == "C1=CC=CC=C1"

def test_sub_search_cache_hit():
    with patch('main.get_cached_result') as mock_get_cached_result:
        mock_get_cached_result.return_value = ["CCO"]

        response = client.get("/subsearch", params={"substructure": "CCO"})
        assert response.status_code == 200
        assert response.json() == {"source": "cache", "data": ["CCO"]}

def test_sub_search_no_cache():
    with patch('main.get_cached_result') as mock_get_cached_result, \
         patch('main.substructure_search_task') as mock_substructure_search_task:

        mock_get_cached_result.return_value = None
        mock_task = MagicMock()
        mock_task.id = "test_task_id"
        mock_substructure_search_task.apply_async.return_value = mock_task

        response = client.get("/subsearch", params={"substructure": "CCO"})
        assert response.status_code == 200
        assert response.json() == {"task_id": "test_task_id", "status": "PENDING"}

def test_get_task_result_pending():
    with patch('main.AsyncResult') as mock_async_result_class:
        mock_async_result_instance = mock_async_result_class.return_value
        mock_async_result_instance.state = 'PENDING'

        response = client.get("/subsearch/task-id")
        assert response.status_code == 200
        assert response.json() == {"task_id": "task-id", "status": "Task is still processing"}

def test_get_task_result_success():
    with patch('main.AsyncResult') as mock_async_result_class:
        mock_async_result_instance = mock_async_result_class.return_value
        mock_async_result_instance.state = 'SUCCESS'
        mock_async_result_instance.result = {'data': 'some result'}

        response = client.get("/subsearch/task-id")
        assert response.status_code == 200
        assert response.json() == {
            "task_id": "task-id",
            "status": "Task completed",
            "result": {'data': 'some result'}
        }



@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    def remove_test_db():
        if os.path.exists("./test.db"):
            os.remove("./test.db")
    request.addfinalizer(remove_test_db)