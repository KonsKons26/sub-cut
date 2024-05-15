from os import listdir
from os.path import isfile, join
from Bio import SeqIO
import pandas as pd
import numpy as np
from itertools import combinations
import multiprocessing as mp


def getData(fileName, path="results/"):
    fastaDict = {}

    processedFiles = [f.rstrip(".csv")
                      for f in listdir(path) if isfile(join(path, f))]

    with open(fileName, "r") as fastaFile:

        for record in SeqIO.parse(fastaFile, "fasta"):
            bigText = record.description
            id, name, meropsID = bigText.strip().split(" | ")
            if "/" in name:
                name = name.replace("/", "-")

            if name not in processedFiles:
                seq = str(record.seq).strip()

                if name not in fastaDict.keys():
                    fastaDict[name] = [seq]

                else:
                    fastaDict[name].append(seq)

    trypsin1 = fastaDict.pop("trypsin 1")
    returnDict = {}
    returnDict["trypsin 1"] = trypsin1

    return returnDict


def blosum62score(x, y):
    blosum62RawTable = [
        [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -
            1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
        [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -
            1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4],
        [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -
            2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4],
        [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -
            3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4],
        [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -
            1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
        [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -
            3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4],
        [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -
            2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
        [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -
            3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4],
        [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -
            2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4],
        [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3,
            1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4],
        [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2,
            2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4],
        [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -
            1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4],
        [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5,
            0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4],
        [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3,
            0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4],
        [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -
            2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
        [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -
            1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4],
        [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -
            1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4],
        [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -
            1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4],
        [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -
            1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4],
        [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2,
            1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4],
        [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -
            3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4],
        [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -
            1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
        [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -
            1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4],
        [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
         -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]

    aminoAcids = ["A", "R", "N", "D", "C", "Q", "E",
                  "G", "H", "I", "L", "K", "M", "F",
                  "P", "S", "T", "W", "Y", "V", "B",
                  "Z", "X", "*"]

    blosum62 = pd.DataFrame(
        blosum62RawTable, columns=aminoAcids, index=aminoAcids)

    return blosum62.loc[x, y]


def staticAligner(seqsList):
    scores = []
    for s1, s2 in combinations(seqsList, 2):
        score = []
        for i in range(len(s1)):
            score.append(blosum62score(s1[i], s2[i]))
        scores.append(score)

    array = np.array(scores)

    means = np.mean(array, axis=0)
    stds = np.std(array, axis=0)

    return means, stds


def processData(name, data):
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(name)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    means, stds = staticAligner(data)
    array = np.vstack((means, stds))
    np.savetxt(f"results/{name}.csv", array, delimiter=",")


def main():
    fastaDict = getData("cleavage_sites.fasta")

    pool = mp.Pool()

    results = [pool.apply_async(processData, args=(k, v))
               for k, v in fastaDict.items()]

    pool.close()
    pool.join()

    for result in results:
        result.get()


if __name__ == "__main__":
    main()
