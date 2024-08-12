# Aligns each combinations of recognition sites for
# each protease as a measure of similarity, using
# biopython's PairwiseAligner

import json
from Bio import Align
from Bio.Align import substitution_matrices
from statistics import mean
from itertools import combinations
import multiprocessing as mp
from collections import OrderedDict

# Load the data in
with open("data/substrates_data.json", "r") as in_file_name:
    substrates_dict = dict(json.load(in_file_name))


def measure_similarity(prot, seqs_list, return_dict):
    """
    Function to measure how similar the substrates are for each protease.
    Each sequence is aligned against each other sequence and their average
    score is calculated and returned in a dictionary. This is the function
    that will be passed to the mp manager.
    """

    print(f"working on: {prot}...")

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("PAM70")
    aligner.mode = "local"
    # Heavily bias our aligner, the less gaps the better when working
    # with protease substrates.
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -5

    scores = []
    # Iterate over each combination of sequences and get their alignment score
    for s1, s2 in combinations(seqs_list, 2):
        aln = aligner.align(s1, s2)
        scores.append(aln.score)

    # The score of proteases with only one substrate is stored as 'None'
    if len(seqs_list) > 1:
        score = mean(scores)
    else:
        score = None

    return_dict[prot] = score


# Multiprocessing magic, calculate the score for each protease in parallel
with mp.Manager() as manager:
    cpus = int(mp.cpu_count())
    # Create the return dictionary
    return_dict = manager.dict()
    # Use one half of the machines cpus
    pool = mp.Pool(processes=int(cpus / 2))
    # Parallel magic
    pool.starmap(measure_similarity, [
                 (prot, seqs_list, return_dict) for prot, seqs_list in substrates_dict.items()])
    pool.close()
    pool.join()

    # OrderdDict so we get the results in alphabetical order
    result = OrderedDict(sorted(return_dict.items()))


others = {}
single_substrate_proteases = []

# Proteases with only one substrate are stored in a different dictionary
for k, v in result.items():
    if v is None:
        single_substrate_proteases.append(k)
    else:
        others[k] = v

# Normalize the values between 0 and 1
max_val = max(others.values())
min_val = min(others.values())
for k, v in others.items():
    # z_i = (x_i – min(x)) / (max(x) – min(x))
    others[k] = (v - min_val) / (max_val - min_val)

# Put both dictionaries in one for saving
final_dict = {}
final_dict["Proteases with more than one substrate"] = others
final_dict["Proteases with only one substrate"] = single_substrate_proteases

# Dump in a json file
with open("data/scored_substrates_biopython.json", "w") as out_file_name:
    json.dump(final_dict, out_file_name, sort_keys=True, indent=4)
