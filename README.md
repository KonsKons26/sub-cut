# sub-cut

Tool to find which proteases cut a substrate based on sequence

## Descreption

The [MEROPS](https://www.ebi.ac.uk/merops/) database is an information resource for peptidases and their substrates. I created this tool to refine the search for proteases that will cut a given peptide. sub-cut will sort all the database's proteases based on the similarity of the provided substrate sequence with the substrates that it's known to cut.

The protease-substrate data from MEROPS db are stored locally but I also provide them processed in the .data/substrates.data.json file.

## Usage

- Pass the sequences that you want cut in `data/target_peptides.fasta`.
- Set the specificity threshold in `sub-cut.py`
- Run `sub-cut.py`

```
python3 sub-cut.py
```

- Check your results in `data/<current_date_and_time>_results.xlsx`
- Check out `analyses.ipynb` for more info and `example.ipynb` for an example
