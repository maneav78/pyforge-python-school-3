import pytest
import os
from databases import Database
from httpx import AsyncClient, ASGITransport
from src.main import app, substructure_search

database = Database(os.getenv("DB_URL"))


@pytest.fixture(scope="module", autouse=True)
async def setup_database():
    await database.connect()
    await database.execute('''
        CREATE TABLE IF NOT EXISTS molecules (
            identifier INTEGER PRIMARY KEY,
            smiles TEXT NOT NULL
        );
    ''')
    yield
    await database.disconnect()


async def get_client():
    return AsyncClient(transport=ASGITransport(app=app), base_url="http://test")


@pytest.mark.asyncio
async def test_substructure_search_valid():
    mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "c1ccccc1"
    result = await substructure_search(mols, substructure)
    assert result == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]


@pytest.mark.asyncio
async def test_substructure_search_invalid():
    mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "invalid_smiles"
    with pytest.raises(ValueError, match='Invalid substructure SMILES: invalid_smiles'):
        await substructure_search(mols, substructure)


async def add_molecule(client, id: int, smiles: str):
    response = await client.post("/add", params={"id": id, "smiles": smiles})
    assert response.status_code == 200
    assert response.json() == {"message": f"Molecule '{id}' added successfully."}


async def get_molecule(client, id: int, expected_smiles: str):
    response = await client.get("/get", params={"id": id})
    assert response.status_code == 200
    assert response.json() == expected_smiles


async def update_molecule(client, id: int, smiles: str):
    response = await client.put("/update", params={"id": id, "smiles": smiles})
    assert response.status_code == 200
    assert response.json() == {"message": f"Molecule '{id}' updated successfully."}


async def delete_molecule(client, id: int):
    response = await client.delete("/del", params={"id": id})
    assert response.status_code == 200
    assert response.json() == {"message": f"Molecule '{id}' deleted successfully."}


async def list_all_molecules(client):
    response = await client.get("/getall")
    assert response.status_code == 200
    return response.json()


async def substructure_search_endpoint(client, substructure: str):
    response = await client.get("/subsearch", params={"substructure": substructure})
    assert response.status_code == 200
    return response.json()


@pytest.mark.asyncio
async def test_add_molecule():
    async with await get_client() as client:
        await add_molecule(client, 1, "CCO")


@pytest.mark.asyncio
async def test_get_molecule():
    async with await get_client() as client:
        await add_molecule(client, 1, "CCO")
        await get_molecule(client, 1, "CCO")


@pytest.mark.asyncio
async def test_update_molecule():
    async with await get_client() as client:
        await add_molecule(client, 1, "CCO")
        await update_molecule(client, 1, "c1ccccc1")
        await get_molecule(client, 1, "c1ccccc1")


@pytest.mark.asyncio
async def test_delete_molecule():
    async with await get_client() as client:
        await add_molecule(client, 1, "CCO")
        await delete_molecule(client, 1)


@pytest.mark.asyncio
async def test_list_all_molecules():
    async with await get_client() as client:
        await add_molecule(client, 1, "CCO")
        await add_molecule(client, 2, "c1ccccc1")
        molecules = await list_all_molecules(client)
        assert len(molecules) == 2


@pytest.mark.asyncio
async def test_substructure_search_endpoint():
    async with await get_client() as client:
        await add_molecule(client, 1, "CCO")
        await add_molecule(client, 2, "c1ccccc1")
        await add_molecule(client, 3, "CC(=O)O")
        matches = await substructure_search_endpoint(client, "c1ccccc1")
        assert "c1ccccc1" in matches
        assert len(matches) == 1
