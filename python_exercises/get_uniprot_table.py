import requests
import pandas as pd

# Get data from the Uniprot using requests


def uniprotquery2pandas(query_text):
    """Send a query to Uniprot and put data in a panda's DataFrame"""
    # query_text = "muscarinic+receptor"
    url = f"http://www.uniprot.org/uniprot/?query={query_text}&sort=score&format=tab"

    req = requests.get(url)
    lines = req.text.splitlines()
    column_names = lines[0].split("\t")
    data_in_rows = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(data=data_in_rows, columns=column_names)
    df.set_index("Entry")
    return df


df = uniprotquery2pandas("opioid")
print("Reviewed median:", df[df["Status"] == "reviewed"]["Length"].median())
print("Uneviewed median:", df[df["Status"] == "unreviewed"]["Length"].median())
