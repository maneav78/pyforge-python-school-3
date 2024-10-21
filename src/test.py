import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from main import app, get_db
from models import Base
from unittest.mock import patch, MagicMock

SQLALCHEMY_TEST_DATABASE_URL = "sqlite:///:memory:"

engine = create_engine(
    SQLALCHEMY_TEST_DATABASE_URL, connect_args={"check_same_thread": False}
)
TestingSessionLocal = sessionmaker(bind=engine)


@pytest.fixture(scope="session")
def db_engine():
    Base.metadata.create_all(bind=engine)
    yield engine
    Base.metadata.drop_all(bind=engine)


@pytest.fixture(scope="function")
def db_session(db_engine):
    connection = db_engine.connect()
    transaction = connection.begin()
    session = TestingSessionLocal(bind=connection)
    app.dependency_overrides[get_db] = lambda: session
    yield session
    session.close()
    transaction.rollback()
    connection.close()


client = TestClient(app)


@pytest.fixture(autouse=True)
def mock_external_dependencies():
    with patch('main.celery'), \
         patch('main.Chem.MolFromSmiles') as mock_mol_from_smiles, \
         patch('main.get_cached_result') as mock_get_cached_result, \
         patch('redis_utils.redis.Redis') as mock_redis_client:

        mock_mol_from_smiles.return_value = MagicMock()

        mock_redis_instance = mock_redis_client.return_value
        mock_redis_instance.get.return_value = None
        yield


@patch('main.Chem.MolFromSmiles')
def test_add_molecule(mock_mol_from_smiles, db_session):
    mock_mol_from_smiles.return_value = MagicMock()

    response = client.post("/add", json={"smiles": "CCO"})
    assert response.status_code == 200
    data = response.json()
    assert data["smiles"] == "CCO"
    assert "id" in data


@patch('main.Chem.MolFromSmiles')
def test_add_molecule_invalid_smiles(mock_mol_from_smiles, db_session):
    mock_mol_from_smiles.return_value = None

    response = client.post("/add", json={"smiles": "InvalidSMILES"})
    assert response.status_code == 400
    data = response.json()
    assert data["detail"] == "Invalid SMILES string"


@patch('main.Chem.MolFromSmiles')
def test_get_molecule(mock_mol_from_smiles, db_session):
    mock_mol_from_smiles.return_value = MagicMock()
    response = client.post("/add", json={"smiles": "C1=CC=CC=C1"})
    data = response.json()
    molecule_id = data["id"]

    response = client.get(f"/get?id={molecule_id}")
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == molecule_id
    assert data["smiles"] == "C1=CC=CC=C1"


@patch('main.Chem.MolFromSmiles')
def test_update_molecule(mock_mol_from_smiles, db_session):
    mock_mol_from_smiles.return_value = MagicMock()
    response = client.post("/add", json={"smiles": "CCO"})
    data = response.json()
    molecule_id = data["id"]

    mock_mol_from_smiles.return_value = MagicMock()

    response = client.put(f"/update?id={molecule_id}", json={"smiles": "C1=CC=CC=C1"})
    assert response.status_code == 200
    data = response.json()
    assert data["smiles"] == "C1=CC=CC=C1"

    response = client.get(f"/get?id={molecule_id}")
    data = response.json()
    assert data["smiles"] == "C1=CC=CC=C1"


@patch('main.Chem.MolFromSmiles')
def test_delete_molecule(mock_mol_from_smiles, db_session):
    mock_mol_from_smiles.return_value = MagicMock()
    response = client.post("/add", json={"smiles": "CCO"})
    data = response.json()
    molecule_id = data["id"]

    response = client.delete(f"/del?id={molecule_id}")
    assert response.status_code == 200
    data = response.json()
    assert data == {"message": f"Molecule '{molecule_id}' deleted successfully."}

    response = client.get(f"/get?id={molecule_id}")
    assert response.status_code == 404


@patch('main.Chem.MolFromSmiles')
def test_list_all_molecules(mock_mol_from_smiles, db_session):
    mock_mol_from_smiles.return_value = MagicMock()
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


@patch('main.get_cached_result', return_value=["CCO"])
def test_sub_search_cache_hit(mock_get_cached_result, db_session):
    response = client.get("/subsearch", params={"substructure": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"source": "cache", "data": ["CCO"]}


@patch('main.get_cached_result', return_value=None)
@patch('main.substructure_search_task')
def test_sub_search_no_cache(mock_substructure_search_task, mock_get_cached_result, db_session):
    mock_task = MagicMock()
    mock_task.id = 'test_task_id'
    mock_substructure_search_task.apply_async.return_value = mock_task

    response = client.get("/subsearch", params={"substructure": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"task_id": "test_task_id", "status": "PENDING"}


@patch('main.AsyncResult')
def test_get_task_result_pending(mock_async_result_class, db_session):
    mock_async_result_instance = mock_async_result_class.return_value
    mock_async_result_instance.state = 'PENDING'

    response = client.get("/subsearch/task-id")
    assert response.status_code == 200
    assert response.json() == {"task_id": "task-id", "status": "Task is still processing"}


@patch('main.AsyncResult')
def test_get_task_result_success(mock_async_result_class, db_session):
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
