# Executes a mysql query from Python (INSERT / DELETE)

import sqlalchemy
import pandas as pd

engine = sqlalchemy.create_engine(
    "mysql+mysqlconnector://fake:omics1819@localhost/fakedb"
)

my_prot = "P12345"
engine.execute(f"INSERT INTO proteins (uniprot_accession) VALUES ('{my_prot}')")
# engine.execute(f"DELETE FROM proteins where uniprot_accession = '{my_prot}'")
