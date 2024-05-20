# The Position Probability Matrix is calculated by counting the
# occurence of each amino acid for each position and the normalizing
# that count between 0 and 1.

import json
import pandas as pd

# Load the data in
with open("data/substrates_data.json", "r") as in_file_name:
    substrates_dict = dict(json.load(in_file_name))


def ppm_creator(sequences_list):
    """
    Create two ppms, one without and one with pseudocounts
    """
    amino_acids = ["*", "X", "A", "R", "N", "D", "C", "Q",
                   "E", "G", "H", "I", "L", "K", "M", "F",
                   "P", "S", "T", "W", "Y", "V"]
    # Create a DataFrame for pfm filled with zeros, where column length = amino acids
    # counts and row length = size of sequences (8)
    pfm = pd.DataFrame({aa: [0] * 8 for aa in amino_acids})
    # Iterate over each sequence
    for sequence in sequences_list:
        # Iterate over each amino acid and add 1 to the corresponding position
        for i, aa in enumerate(sequence):
            pfm.loc[i, aa] += 1
    # Normalize to get the ppm
    ppm = pfm.apply(lambda row: row / row.sum(), axis=1)

    return ppm


# Loop over all proteases and add the results to dictionaries
ppm_dict = {}
for k, v in substrates_dict.items():
    ppm = ppm_creator(v)
    ppm_dict[k] = ppm.to_dict(orient="list")

# Dump in a json file
with open("data/ppm.json", "w") as out_file_name:
    json.dump(ppm_dict, out_file_name, sort_keys=True, indent=4)
