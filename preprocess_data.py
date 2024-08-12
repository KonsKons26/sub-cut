# The merops DB is stored locally in MySQL, the
# required data is fetched and sorted out,
# then saved in json format.

import configparser
from sqlalchemy import create_engine
import pandas as pd
import json


# Personal details hidden in a .config.ini file which is included in .gitignore
config = configparser.ConfigParser()
config.read(".config.ini")

# Set up the SQLAlchemy engine
DB_USER = config.get("database", "user")
DB_PASSWORD = config.get("database", "password")
DB_HOST = config.get("database", "host")
DB_PORT = config.get("database", "port")
DB_NAME = config.get("database", "database")
db_uri = f"mysql+pymysql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
engine = create_engine(db_uri)

# SQL query
query = """
    SELECT
        Uniprot,
        code,
        Protease,
        Site_P4,
        Site_P3,
        Site_P2,
        Site_P1,
        Site_P1prime,
        Site_P2prime,
        Site_P3prime,
        Site_P4prime
    FROM Substrate_search
    WHERE
        Uniprot IS NOT NULL
        AND code IS NOT NULL
"""

# Store data in a DF
df = pd.read_sql(query, engine)

# Empty column that will be filled with the formatted sequences
df["Substrate_Sequence"] = None


def convert_aa(aa):
    """
    Function to convert the three letter code amino acids to single letter.
    As seen from Bio/Align/substitution_matrices/data/<matrix>,
    '*' is used as a gap character (for PAM and BLOSUM atleast).
    """
    aa_dict = {"xaa": "X", "ala": "A", "arg": "R", "asn": "N", "asp": "D", "cys": "C", "gln": "Q",
               "glu": "E", "gly": "G", "his": "H", "ile": "I", "leu": "L", "lys": "K", "met": "M",
               "phe": "F", "pro": "P", "ser": "S", "thr": "T", "trp": "W", "tyr": "Y", "val": "V"}
    if aa.lower() in aa_dict.keys():
        new_aa = aa_dict[aa.lower()]
    else:
        new_aa = '*'

    return new_aa


# Iterate over each row and format the substrate sequence, then add it to the empty column
for i, col in df.iterrows():
    # 3 letter code
    letters3 = col[3:-1]
    # 1 letter code
    letters1 = "".join([convert_aa(letter) for letter in letters3])
    df.loc[i, "Substrate_Sequence"] = letters1

# Group substrate sequences by protease name and place in a dictionary
d = df.groupby("Protease")["Substrate_Sequence"].apply(list).to_dict()

# Convert the lists of sequences to sets to remove any
# duplicates, back again to lists so they work with json,
# and overwrite them in the dictionary
bigger_than_5k = []
new_d = {}
for k, v in d.items():
    s = list(set(v))
    if len(s) < 5_000:
        new_d[k] = s
    else:
        bigger_than_5k.append(k)

print(
    f"Proteases {bigger_than_5k} were ignored, as\nthey have more than 5k substrate sequences.")

# Dump in a json file
with open("data/substrates_data.json", "w") as out_file_name:
    json.dump(new_d, out_file_name, sort_keys=True, indent=4)
