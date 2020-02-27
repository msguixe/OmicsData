# Executes a mysql query from Python and brings the results to a pandas dataframe

import sqlalchemy
import pandas as pd

engine = sqlalchemy.create_engine(
    "mysql+mysqlconnector://fake:omics1819@localhost/fakedb"
)

df = pd.read_sql(
    "SELECT uniprot_entry, gene_name FROM proteins", engine, index_col="uniprot_entry"
)
