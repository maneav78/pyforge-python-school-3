import pytest
from src.main import substructure_search, app, create_table
from fastapi.testclient import TestClient
import time
import os

client = TestClient(app)


@pytest.fixture(scope="module", autouse=True)
def setup_database():
    if os.path.exists("molecules.db"):
        os.remove("molecules.db")
    create_table()
    yield
    if os.path.exists("molecules.db"):
        os.remove("molecules.db")


def test_substructure_search_valid():
    mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "c1ccccc1"
    result = substructure_search(mols, substructure)
    assert result == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]


def test_substructure_search_invalid():
    mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "invalid_smiles"
    with pytest.raises(ValueError, match='Invalid substructure SMILES: invalid_smiles'):
        substructure_search(mols, substructure)


def test_add_molecule():
    time.sleep(2)
    response = client.post("/add", params={"id": 1, "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule '1' added successfully."}


def test_get_molecule():
    response = client.post("/add", params={"id": 2, "smiles": "c1ccccc1"})
    assert response.status_code == 200
    response = client.get("/get", params={"id": 2})
    assert response.status_code == 200
    assert response.json() == "c1ccccc1"


def test_update_molecule():
    response = client.post("/add", params={"id": 3, "smiles": "CC(=O)O"})
    assert response.status_code == 200
    response = client.put("/update", params={"id": 3, "smiles": "c1ccccc1"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule '3' updated successfully."}
    response = client.get("/get", params={"id": 3})
    assert response.status_code == 200
    assert response.json() == "c1ccccc1"


def test_delete_molecule():
    response = client.post("/add", params={"id": 4, "smiles": "CC(=O)Oc1ccccc1C(=O)O"})
    assert response.status_code == 200
    response = client.delete("/del", params={"id": 4})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule '4' deleted successfully."}
    response = client.get("/get", params={"id": 4})
    assert response.status_code == 404


def test_list_all_molecules():
    response = client.post("/add", params={"id": 5, "smiles": "CCO"})
    assert response.status_code == 200
    response = client.post("/add", params={"id": 6, "smiles": "c1ccccc1"})
    assert response.status_code == 200
    response = client.get("/getall")
    assert response.status_code == 200


def test_substructure_search_endpoint():
    response = client.post("/add", params={"id": 7, "smiles": "CCO"})
    assert response.status_code == 200
    response = client.post("/add", params={"id": 8, "smiles": "c1ccccc1"})
    assert response.status_code == 200
    response = client.get("/subsearch", params={"substructure": "c1ccccc1"})
    assert response.status_code == 200
    assert "c1ccccc1" in response.json()
