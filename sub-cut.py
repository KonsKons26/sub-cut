import json
import pandas as pd
import multiprocessing as mp
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from datetime import datetime


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

    return aligner


def get_params(score_type="positionwise", include_proteases_with_one_substrate=True):
    # Load the data
    with open("data/substrates_data.json", "r") as in_file_name:
        substrates_dict = dict(json.load(in_file_name))

    if score_type == "positionwise":
        with open("data/scored_substrates_positionwise.json", "r") as scores_file:
            scores_dict = dict(json.load(scores_file))
    elif score_type == "biopython":
        with open("data/scored_substrates_biopython.json", "r") as scores_file:
            scores_dict = dict(json.load(scores_file))
    else:
        raise ValueError(
            "Invalid score_type: must be 'positionwise' or 'biopython'")

    # Set the specificity threshold here:
    # range from 0 to 1, 0 will include all protease,
    # 0.25 will include the top 75% of protease,
    # 0.42 will include the top 58% of proteases,
    # 0.5 will include the top 50% of proteases,
    # 0.8 will include the top 20% of proteases, etc
    thresh = 0.6

    proteases_dict = {}

    # Always include proteases with only one substrate, as they are highly specific,
    # unless otherwise specified
    for k, v in substrates_dict.items():
        if include_proteases_with_one_substrate:
            for one_substrate_proteases in scores_dict["Proteases with only one substrate"]:
                proteases_dict[one_substrate_proteases] = substrates_dict[one_substrate_proteases]
        # exclude the cursed trypsin 1
        if k not in scores_dict["Proteases with only one substrate"] and k != "trypsin 1":
            if scores_dict["Proteases with more than one substrate"][k] >= thresh:
                proteases_dict[k] = v

    targets = {}
    #  Get the target name and target_sequence from the .txt file
    with open("data/target_peptides.fasta", "r") as targets_file:
        records = SeqIO.parse(targets_file, "fasta")
        for record in records:
            targets[str(record.id)] = str(record.seq)

    return proteases_dict, targets


def batchalign(proteases_dict, aligner, target_name, target_seq):
    # Create empty lists
    proteases = []
    rec_sites = []
    scores = []
    alignments = []

    # Loop over each protease, substrates_list pair
    for k, v in proteases_dict.items():
        # Loop over each substrate sequence
        for seq in v:
            # This will look like this:
            # ["Prot1", "Prot1", "Prot1", "Prot2", "Prot3",...]
            # ["P1seq1", "P1seq2", "P1seq3", "P2seq1", "P3seq1",...]
            proteases.append(k)
            rec_sites.append(seq)
            try:
                # Align the protease substrate seq with the target peptide seq
                aln = aligner.align(target_seq, seq)
                score = aln.score
                alignment = str(aln[0]).rstrip("\n")
            except IndexError:
                # IndexError raised on `aln[0]` if no alignments are made
                score = None
                alignment = "No significant alignment was made!"
            scores.append(score)
            alignments.append(alignment)

    # Store all lists in a DF
    df = pd.DataFrame(
        {
            "Protease": proteases,
            "Recognition Site": rec_sites,
            "Score": scores,
            "Alignment": alignments}
    )

    #  Sort the scores by descending order
    df = df.sort_values(by='Score', ascending=False)

    # Headers must have the form:
    # >7BZ5_1|Chain_A_433-510|Spike_protein_S1
    # for this to work properly, the file name will look like this:
    # <currentDate&Time>-7BZ5_1Chain_A_433-510.xlsx
    peptide_code = "".join(target_name.split("|")[:2])
    now = datetime.now()
    out_file_name = f"results/{now.year}-{now.month}-{now.day}_{now.hour}{now.minute}-{peptide_code}.xlsx"

    # Write the output in an xlsx file so that the alignments can also be visualized
    with pd.ExcelWriter(out_file_name, engine="xlsxwriter") as xlsxwriter:
        df.to_excel(xlsxwriter, sheet_name="Sheet1", index=False)
        workbook = xlsxwriter.book
        format = workbook.add_format(
            {
                "text_wrap": True,
                "font_name": "Consolas"
            }
        )
        worksheet = xlsxwriter.sheets["Sheet1"]
        worksheet.set_column("A:F", None, cell_format=format)


def multiprocessing(proteases_dict, aligner, targets_dict):
    # MP magic, run every target peptide in parallel
    processes = []
    for target_name, target_seq in targets_dict.items():

        p = mp.Process(target=batchalign, args=[
            proteases_dict,
            aligner,
            target_name,
            target_seq
        ])
        p.start()
        processes.append(p)

    [process.join() for process in processes]


def move_targets():
    # Appends the target peptides processed to the file data/old_target_peptides.txt
    # and empties out the data/target_peptides.fasta
    now = datetime.now()
    now = f"{now.year}-{now.month}-{now.day}_{now.hour}{now.minute}"
    with open("data/target_peptides.fasta", "r") as current_targets_file:
        current_targets = current_targets_file.read()
    with open("data/old_target_peptides.txt", "a") as old_targets_file:
        old_targets_file.write(f"{now}\n")
        old_targets_file.write(current_targets)
        old_targets_file.write("\n")
    with open("data/target_peptides.fasta", "w") as current_targets_file:
        current_targets_file.write("")


def main():
    #  Run the functions
    aligner = setup_aligner()
    proteases_dict, targets_dict = get_params()
    multiprocessing(proteases_dict, aligner, targets_dict)
    move_targets()


if __name__ == '__main__':
    main()
