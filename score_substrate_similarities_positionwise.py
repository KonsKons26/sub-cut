# Aligns each combination of recognition sites for
# each protease using a custom static alignment
# method that checks each position with each corresponding
# position(x1 - y1, x2 - y2, ..., x8 - y8) to measure
# the degree of conservation of each site.
# Scoring is based on the BLOSUM62 table.

import json
import pandas as pd
from itertools import combinations
import multiprocessing as mp

# Load the data in
with open("data/substrates_data.json", "r") as in_file_name:
    substrates_dict = dict(json.load(in_file_name))

# Initialize new dicts
one_substrate_proteases = []
mul_substrate_proteases = {}

# Separate proteases with only one substrate
for k, v in substrates_dict.items():
    if len(v) == 1:
        one_substrate_proteases.append(k)
    else:
        mul_substrate_proteases[k] = v

# Set cpus to be used, we will split the data in a list of int(mp.cpu_count() / 1.25) smaller dicts
cpus = int(mp.cpu_count() / 1.25)
# Initialize the list of smaller dicts
smaller_dicts = [{} for _ in range(cpus)]
# Keep track of current required operations for each index of the smaller_dicts list
current_ops_count = [0] * cpus

# Loop over k: v pairs of proteases
for k, v in mul_substrate_proteases.items():
    # Calculate required operations for this protease
    op_count = int((len(v) * (len(v) - 1)) / 2)
    # Check which index has the least operations assigned
    min_ops_index = current_ops_count.index(min(current_ops_count))
    # Add the smaller_dict to the index with the least operations
    smaller_dicts[min_ops_index][k] = v
    # Update the list of operations for each index
    current_ops_count[min_ops_index] += op_count

# Print some general data to quickly check that everythin is running smoothly
for i, d in enumerate(smaller_dicts):
    print(
        f"Smaller Dictionary {i+1}\tn. of operations:\t{sum([((len(v) * len(v) - 1) / 2) for v in d.values()])}")
print("Missing data:",
      len(substrates_dict) - sum([len(d.values())
                                  for d in smaller_dicts]) - len(one_substrate_proteases))

# Blosum62 table data
blosum62df = pd.read_csv("data/blosum62.csv",  index_col=0)


def measure_similarity(return_list, blosum62, smaller_dict):
    """
    Function to measure how similar the substrates are for each protease.
    Each sequence is aligned against each other sequence and their average
    score is calculated and returned in a dictionary. This is the function
    that will be passed to the mp manager.
    """
    # Initialize empty list
    l = []
    # Loop over each protease - substrate sequences list pair
    # and calculate its score positionwise
    for protease, seqs_list in smaller_dict.items():
        print(f"working on: {protease}...")
        protease_scores = []
        for seq1, seq2 in combinations(seqs_list, 2):
            position_scores = []
            for aa1, aa2 in zip(seq1, seq2):
                position_scores.append(blosum62.loc[aa1, aa2])
            protease_scores.append(sum(position_scores) / 8)
        score = sum(protease_scores) / len(protease_scores)
        l.append({protease: score})

    return_list.append(l)


# Multiprocessing magic, calculate the score for each protease in parallel
with mp.Manager() as manager:
    # Initialize a return list for appending the resulting dicts
    return_list = manager.list()
    # Start processing in parallel using the same number of cpus calculated above
    pool = mp.Pool(processes=cpus)
    pool.starmap(
        measure_similarity, [
            (return_list, blosum62df, smaller_dict) for smaller_dict in smaller_dicts
        ]
    )
    pool.close()
    pool.join()
    # Conver the resulting list to a proper list, before closing the Manager()
    result = list(return_list)


def flatten_dict_list(nested_list):
    """
    Recursive function to flatten the output list and make it a dict
    of {protease: score} pairs
    """
    flattened_dict = {}
    for item in nested_list:
        # Loop over each item in the nested list
        if isinstance(item, dict):
            # If its a dict, add it to the final flattened dict
            flattened_dict.update(item)
        elif isinstance(item, list):
            # If its a list, call the function again
            flattened_dict.update(flatten_dict_list(item))
    return flattened_dict


# Flatten the nested list and convert to a single dictionary, then order
result_dict = flatten_dict_list(result)
result_dict = dict(sorted(result_dict.items(), key=lambda x: x[1]))

# Normalize the values between 0 and 1
max_val = max(result_dict.values())
min_val = min(result_dict.values())
max_minus_min = max(result_dict.values()) - min(result_dict.values())

for k, v in result_dict.items():
    # z_i = (x_i – min(x)) / (max(x) – min(x))
    result_dict[k] = (v - min_val) / (max_minus_min)

# Put both dictionaries in one for saving
final_dict = {}
final_dict["Proteases with more than one substrate"] = result_dict
final_dict["Proteases with only one substrate"] = one_substrate_proteases

# Dump in a json file
with open("data/scored_substrates_positionwise.json", "w") as out_file_name:
    json.dump(final_dict, out_file_name, sort_keys=True, indent=4)
