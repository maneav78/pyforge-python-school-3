from fastapi import FastAPI, Depends, HTTPException, File, UploadFile, Query
from rdkit import Chem
import sqlite3
from sqlite3 import Error
from os import getenv


app = FastAPI()

db = "molecules.db"
# DEFAULT_FILENAME = "uploaded_image.pdf"


def substructure_search(mols, mol):
    substructure = Chem.MolFromSmiles(mol)

    if substructure is None:
        raise ValueError(f'Invalid substructure SMILES: {mol}')

    substructure_num_atoms = substructure.GetNumAtoms()

    matches = []
    for molecule in mols:
        object_mol = Chem.MolFromSmiles(molecule)
        if object_mol is None:
            raise ValueError(f'Invalid Molecule!')

        if object_mol.GetNumAtoms() < substructure_num_atoms:
            continue

        if object_mol.HasSubstructMatch(substructure):
            matches.append(molecule)

    return matches


def get_connection():
    connection = None
    try:
        connection = sqlite3.connect(db)
        print(f"Connected to SQLite database: {sqlite3.version}")
        return connection
    except Error as e:
        print(e)
    return connection


def create_table():
    connection = get_connection()
    sql_code = """
    CREATE TABLE IF NOT EXISTS molecules (
        identifier INTEGER PRIMARY KEY,
        smiles TEXT NOT NULL
    )
    """
    try:
        c = connection.cursor()
        c.execute(sql_code)
        connection.commit()
        print("Table 'molecules' created successfully.")
    except Error as e:
        print(e)
    finally:
        if connection:
            connection.close()


@app.on_event("startup")
def startup_event():
    create_table()

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

@app.post("/add")
def add_molecule(id: int, smiles: str, connection=Depends(get_connection)):
    sql_code = """
    INSERT INTO molecules (identifier, smiles)
    VALUES (?, ?)
    """
    try:
        c = connection.cursor()
        c.execute(sql_code, (id, smiles))
        connection.commit()
        return {"message": f"Molecule '{id}' added successfully."}
    except Error as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/get")
def get_molecule(id: int, connection=Depends(get_connection)):
    sql_code = """
    SELECT smiles FROM molecules WHERE identifier = ?
    """
    try:
        c = connection.cursor()
        c.execute(sql_code, (id,))
        row = c.fetchone()
        if row:
            return row[0]
        else:
            raise HTTPException(status_code=404, detail="Molecule not found")
    except Error as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.put("/update")
def update_molecule(id: int, smiles: str, connection=Depends(get_connection)):
    sql_code = """
    UPDATE molecules SET smiles = ? WHERE identifier = ?
    """
    try:
        c = connection.cursor()
        c.execute(sql_code, (smiles, id))
        connection.commit()
        return {"message": f"Molecule '{id}' updated successfully."}
    except Error as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.delete("/del")
def delete_molecule(id: int, connection=Depends(get_connection)):
    sql_code = """
    DELETE FROM molecules WHERE identifier = ?
    """
    try:
        c = connection.cursor()
        c.execute(sql_code, (id,))
        connection.commit()
        return {"message": f"Molecule '{id}' deleted successfully."}
    except Error as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/getall")
def list_all_molecules(connection=Depends(get_connection)):
    sql_code = """
    SELECT identifier, smiles FROM molecules
    """
    try:
        c = connection.cursor()
        c.execute(sql_code)
        rows = c.fetchall()
        return rows
    except Error as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/subsearch")
def sub_search(substructure: str, connection=Depends(get_connection)):
    sql_code = """
    SELECT smiles FROM molecules
    """
    try:
        c = connection.cursor()
        c.execute(sql_code)
        rows = c.fetchall()
        molecules = [row[0] for row in rows]
        matches = substructure_search(molecules, substructure)
        return matches
    except Error as e:
        raise HTTPException(status_code=500, detail=str(e))


# Upload file comming soon...

# Testing code
if __name__ == "__main__":
    # Test with valid substructure
    assert substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1") == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]

    # Test with invalid substructure
    try:
        substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "invalid_smiles")
    except ValueError as e:
        assert str(e) == 'Invalid substructure SMILES: invalid_smiles'

    print("All tests passed!")
