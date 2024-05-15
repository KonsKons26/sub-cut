
def setup_aligner():

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

    return aligner, matrix


def get_params():

    thresh = 17.5
    file_name = f"filtered_data_{thresh}.csv"

    #  Get the target name and sequence from the .txt file
    targets = {}
    with open("targets_dictionary.txt") as f:
        for line in f:
            (key, val) = line.split(':')
            val = val.rstrip()
            targets[key] = val

        my_path = os.path.dirname(os.path.abspath(__file__))

    return file_name, targets, my_path, thresh


def batchalign(file_name, target, aligner, my_path, matrix, targetname, thresh):

    #  Read the formatted file
    df = pd.read_csv(file_name)

    #  Save the recognition sites and peptide names to different lists
    proteases = df['Protease'].tolist()
    sites1 = df['recognition_site'].tolist()
    uniprot = df['UniProt'].tolist()

    sites = []

    #  Convert the sequences to match the biopython alphabet
    for site in sites1:
        nsite = ''
        for letter in site:
            if letter == '-':
                nsite += 'X'
            else:
                nsite += letter
        sites.append(nsite)

    #  Create empty lists for scores
    scores = []
    alignments = []
    #  Iterate over every recognition site and perform the alignment
    #  Then append the results to the lists
    for seq in sites:
        try:
            aln = aligner.align(target, seq)
            scores.append(aln.score)
            alignments.append(str(aln[0]).rstrip('\n'))
        except IndexError:
            alignments.append('No significant alignment was made!')

    #  Create a new Data Frame to save all the data
    new_df = pd.DataFrame(
        {
            'Protease': proteases,
            'UniProt': uniprot,
            'Rec_site': sites,
            'Score': scores,
            'Alignment': alignments
        })

    #  Sort the scores by descending order
    final_df = new_df.sort_values(by='Score', ascending=False)

    #  Check that a directory exists for the currently used martix
    Path(my_path + '/results/thresh{}/{}'.format(thresh, matrix)).mkdir(
        parents=True, exist_ok=True)

    #  Save the results
    with pd.ExcelWriter(
            my_path +
        '/results/thresh{}/{}/{}.xlsx'.format(thresh,
                                              matrix, targetname, sep=','),
            engine='xlsxwriter') as writer:

        final_df.to_excel(writer, sheet_name='Sheet1', index=False)
        workbook = writer.book
        format = workbook.add_format(
            {'text_wrap': True,
             'font_name': 'Consolas'})
        worksheet = writer.sheets['Sheet1']
        worksheet.set_column('A:F', None, cell_format=format)

    return


def multiprocessing(file_name, targets, my_path, aligner, matrix, thresh):

    processes = []
    for targetname, target in targets.items():

        p = mp.Process(target=batchalign, args=[
            file_name,
            target,
            aligner,
            my_path,
            matrix,
            targetname,
            thresh
        ])
        p.start()
        processes.append(p)

    [process.join() for process in processes]

    return


def main():

    #  Run the functions
    aligner, matrix = setup_aligner()
    file_name, targets, my_path, thresh = get_params()
    multiprocessing(file_name, targets, my_path, aligner, matrix, thresh)

    return


if __name__ == '__main__':

    import os
    import pandas as pd
    import multiprocessing as mp
    from Bio import Align
    from Bio.Align import substitution_matrices
    from pathlib import Path

    usage = """
    Reads the filtered data and aligns each sequence
    from the filtered recognition sites with each
    sequence from the 'targets_dictionary.txt' file.
    Saves the results in an .xlsx file in a new directory.
    """
    print(usage)

    main()
