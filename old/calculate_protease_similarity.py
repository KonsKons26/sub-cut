
def load_data():

    with open('sequences.p', 'rb') as pf:
        seqs = pickle.load(pf)

    return seqs


def measure_similarity(prot, seqs_list, return_dict):

    #  Set up the aligner
    aligner = Align.PairwiseAligner()

    # Valid Substitution Matrices
    # 'BENNER22', 'BENNER6', 'BENNER74', 'BLOSUM45',
    # 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90',
    # 'DAYHOFF', 'FENG', 'GENETIC', 'GONNET1992', 'HOXD70',
    # 'JOHNSON', 'JONES', 'LEVIN', 'MCLACHLAN', 'MDM78',
    # 'NUC.4.4', 'PAM250', 'PAM30', 'PAM70', 'RAO',
    # 'RISLER', 'SCHNEIDER', 'STR', 'TRANS'
    matrix = 'PAM70'

    aligner.substitution_matrix = substitution_matrices.load(matrix)
    aligner.mode = 'local'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -5

    sim = 0
    counter = 0
    for s1, s2 in combinations(seqs_list, 2):
        aln = aligner.align(s1, s2)
        sim += aln.score
        counter += 1

    if counter == 0:
        avg_sim = 60
    else:
        avg_sim = sim / counter

    return_dict[prot] = avg_sim


def main():

    seqs = load_data()

    with mp.Manager() as manager:
        cpus = int(mp.cpu_count())
        return_dict = manager.dict()
        pool = mp.Pool(processes=int(cpus / 2))
        pool.starmap(measure_similarity, [
                     (prot, seqs_list, return_dict) for prot, seqs_list in seqs.items()])
        pool.close()
        pool.join()

        result = OrderedDict(sorted(return_dict.items()))
        result = list(result.items())

        df = pd.DataFrame(
            result, columns=['Protease', 'Normalized similarity based on alignments'])
        df.to_csv('normalized_ali_score.csv', index=False)


if __name__ == '__main__':
    import pickle
    import pandas as pd
    from itertools import combinations
    from Bio import Align
    from Bio.Align import substitution_matrices
    import multiprocessing as mp
    from collections import OrderedDict

    usage = """
    Aligns each combinations of recognition sites for
    each protease as a measure of similarity and saves
    the resulting file in 'normalized_ali_score.csv'.
    """
    print(usage)

    main()
