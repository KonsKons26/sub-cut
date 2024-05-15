# The Position Probability Matrix is calculated by counting the
# occurence of each amino acid for each position and the normalizing
# that count between 0 and 1.

# The Position Specific Scoring Matrix (PSSM) is calculated as follows:
# M_ik = log_2(M_ik/b_k), where the element at position ik (M_ik)
# equals the log_2 of the element at ik divided by the likelyhood of
# it being k. b_k = 1/|k| (so 0.25 for nucleotides and 0.05 for
# amino acids). Pseudocounts will be added to avoid negative infinities.
# When the PSSM elements are calculated using log likelihoods, the score
# of a sequence can be calculated by adding (rather than multiplying) the
# relevant values at each position in the PSSM. The sequence score gives an
# indication of how different the sequence is from a random sequence. The
# score is 0 if the sequence has the same probability of being a functional
# site and of being a random site. The score is greater than 0 if it is more
# likely to be a functional site than a random site, and less than 0 if
# it is more likely to be a random site than a functional site.

# Reads the substrate sequences and:
# 1. calculates a position probability matrix and stores the data
# in a json file as a dictionary of dictionaries of lists
# {"Protease": {"amino acid": [position probabilites], ...}, ...}

# 2. calculates 3 pssm's. The first pssm will have equal likelyhoods for
# every amino acid, the second one will have amino acid likelyhoods calculated
# separately for each protease and the third one will have likelyhoods
# calculated from the whole dataset.


import json
import pandas as pd
from math import log2


# Load the data in
with open("data/substrates_data.json", "r") as in_file_name:
    substrates_dict = dict(json.load(in_file_name))


def all_likelyhoods(d):
    """
    Returns three dictionaries of letter likelyhoods.
    One with static likelyhoods of the form:
    {letter: likelyhood, ...}
    One calculated for each protease of the form:
    {protease: {letter: likelyhood, ...}, ...}
    One calculated using the whole dataset of the from:
    {letter: likelyhood, ...}
    """
    likelyhoods = {"*": 1 / 22, "X": 1 / 22, "A": 1 / 22, "R": 1 / 22, "N": 1 / 22,
                   "D": 1 / 22, "C": 1 / 22, "Q": 1 / 22, "E": 1 / 22, "G": 1 / 22,
                   "H": 1 / 22, "I": 1 / 22, "L": 1 / 22, "K": 1 / 22, "M": 1 / 22,
                   "F": 1 / 22, "P": 1 / 22, "S": 1 / 22, "T": 1 / 22, "W": 1 / 22,
                   "Y": 1 / 22, "V": 1 / 22}
    # N of all letters (aas and gaps) for the whole dataset calculation of likelyhoods
    letters_sum_dataset = sum(len(inner_list) * 8 for inner_list in d.values())
    # Use pseudocounts equal to 1% of dataset size
    pseudocounts = .01
    # Create the empty map for the dataset calculated likelyhoods
    likelyhoods_dataset = {"*": pseudocounts, "X": pseudocounts, "A": pseudocounts, "R": pseudocounts,
                           "N": pseudocounts, "D": pseudocounts, "C": pseudocounts, "Q": pseudocounts,
                           "E": pseudocounts, "G": pseudocounts, "H": pseudocounts, "I": pseudocounts,
                           "L": pseudocounts, "K": pseudocounts, "M": pseudocounts, "F": pseudocounts,
                           "P": pseudocounts, "S": pseudocounts, "T": pseudocounts, "W": pseudocounts,
                           "Y": pseudocounts, "V": pseudocounts}
    likelyhoods_per_protease = {}
    # Loop over each protease, sequences list pair
    for protease, seqs_list in d.items():
        # N of all letters (aas and gaps) for the per protease calculation of likelyhoods
        letters_sum_protease = 8 * len(seqs_list)
        # Use pseudocounts equal to 1% of dataset size
        pseudocounts = .01
        # Create the empty map for the likelyhoods per protease
        likelyhoods_per_protease[protease] = {"*": pseudocounts, "X": pseudocounts, "A": pseudocounts, "R": pseudocounts,
                                              "N": pseudocounts, "D": pseudocounts, "C": pseudocounts, "Q": pseudocounts,
                                              "E": pseudocounts, "G": pseudocounts, "H": pseudocounts, "I": pseudocounts,
                                              "L": pseudocounts, "K": pseudocounts, "M": pseudocounts, "F": pseudocounts,
                                              "P": pseudocounts, "S": pseudocounts, "T": pseudocounts, "W": pseudocounts,
                                              "Y": pseudocounts, "V": pseudocounts}
        # Loop over each sequence
        for seq in seqs_list:
            # Loop over each aa
            for aa in seq:
                # Add one to the corresponding place
                # in the maps and to the total letters sums
                likelyhoods_dataset[aa] += 1
                letters_sum_dataset += 1
                likelyhoods_per_protease[protease][aa] += 1
                letters_sum_protease += 1
        # Loop over the protease map and divide the counts vy the total
        # to get the likelyhood of each letter
        likelyhoods_per_protease[protease] = {letter1: n1 / letters_sum_protease for letter1, n1
                                              in likelyhoods_per_protease[protease].items()}
    # Loop over the whole dataset map and divide the counts by the total
    # to get the likelyhood of each letter
    likelyhoods_dataset = {letter2: n2 / letters_sum_dataset for letter2, n2
                           in likelyhoods_dataset.items()}

    return likelyhoods, likelyhoods_per_protease, likelyhoods_dataset


def ppm_creator(sequences_list):
    """
    Create two ppms, one without and one with pseudocounts
    """
    amino_acids = ["*", "X", "A", "R", "N", "D", "C", "Q",
                   "E", "G", "H", "I", "L", "K", "M", "F",
                   "P", "S", "T", "W", "Y", "V"]
    # Create a DataFrame for pfm filled with zeros, where column length = amino acids
    # counts and row length = size of sequences (8) and a pfm with pseudocounts
    pseudocounts = .01
    pfm = pd.DataFrame({aa: [0] * 8 for aa in amino_acids})
    pfm_with_pseudocounts = pd.DataFrame(
        {aa: [pseudocounts] * 8 for aa in amino_acids})
    # Iterate over each sequence
    for sequence in sequences_list:
        # Iterate over each amino acid and add 1 to the corresponding position
        for i, aa in enumerate(sequence):
            pfm.loc[i, aa] += 1
            pfm_with_pseudocounts.loc[i, aa] += 1
    # Normalize to get the ppm
    ppm = pfm.apply(lambda row: row / row.sum(), axis=1)
    ppm_with_pseudocounts = pfm_with_pseudocounts.apply(
        lambda row: row / row.sum(), axis=1)

    return ppm, ppm_with_pseudocounts


def pssm_creator(ppm_pseudocounts_dict, likelyhoods, likelyhoods_per_protease, likelyhoods_dataset):
    """
    Create the three pssm matrices
    """
    pssm = {}
    pssm_protease = {}
    pssm_dataset = {}
    for protease, df in ppm_pseudocounts_dict.items():
        pssm[protease] = df.apply(lambda col:
                                  col.apply(lambda x:
                                            log2(x / likelyhoods[col.name]))
                                  ).to_dict(orient="list")
        pssm_protease[protease] = df.apply(lambda col:
                                           col.apply(lambda x:
                                                     log2(x / likelyhoods_per_protease[protease][col.name]))
                                           ).to_dict(orient="list")
        pssm_dataset[protease] = df.apply(lambda col:
                                          col.apply(lambda x:
                                                    log2(x / likelyhoods_dataset[col.name]))
                                          ).to_dict(orient="list")

    return pssm, pssm_protease, pssm_dataset


# Calculate the letter likelyhoods
likelyhoods, likelyhoods_per_protease, likelyhoods_dataset = all_likelyhoods(
    substrates_dict)

# Loop over all proteases and add the results to dictionaries
ppm_dict = {}
ppm_pseudocounts_dict = {}
for k, v in substrates_dict.items():
    ppm, ppm_pseudocounts = ppm_creator(v)
    ppm_dict[k] = ppm.to_dict(orient="list")
    ppm_pseudocounts_dict[k] = ppm_pseudocounts

# Calculate the three pssm's
pssm_dict, pssm_protease_dict, pssm_dataset_dict = pssm_creator(ppm_pseudocounts_dict, likelyhoods,
                                                                likelyhoods_per_protease, likelyhoods_dataset)

# Dump in a json file
with open("data/ppm.json", "w") as out_file_name:
    json.dump(ppm_dict, out_file_name, sort_keys=True, indent=4)

with open("data/pssm.json", "w") as pssm1:
    json.dump(pssm_dict, pssm1, sort_keys=True, indent=4)

with open("data/pssm_protease.json", "w") as pssm2:
    json.dump(pssm_protease_dict, pssm2, sort_keys=True, indent=4)

with open("data/pssm_dataset.json", "w") as pssm3:
    json.dump(pssm_dataset_dict, pssm3, sort_keys=True, indent=4)
