def load_data():

    data1 = pd.read_csv('cleavage_sites_UniProt_clean.csv')

    data2 = pd.read_csv('normalized_ali_score.csv', header=0, index_col=0)
    data2['Normalized similarity based on alignments'] = data2['Normalized similarity based on alignments'].astype(
        float)

    return data1, data2


def plot_diff_dist(data):

    d = data['Normalized similarity based on alignments'].tolist()

    _, ax = plt.subplots()
    ax.hist(d, bins=750)
    ax.set_yscale('log')
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_title('Normalized similarity score based on alignments',
                 fontdict={
                     'fontsize': 18,
                     'fontweight': 7.5,
                 })
    ax.set_xlabel('Normalized similarity',
                  fontdict={
                      'fontsize': 14,
                      'fontweight': 5,
                  })
    ax.set_ylabel('Counts',
                  fontdict={
                      'fontsize': 14,
                      'fontweight': 5,
                  })
    plt.show()
    return


def apply_mask(data, mask):

    thresh = 17.5
    clean_data = pd.merge(data, mask, how='left', on='Protease')
    clean_data = clean_data.loc[clean_data['Normalized similarity based on alignments'] > thresh]

    clean_data.to_csv(f"filtered_data_{thresh}.csv")

    return clean_data


def main():

    data1, data2 = load_data()
    # plot_diff_dist(data2)  # Uncomment to create the plot
    apply_mask(data1, data2)

    return


if __name__ == '__main__':
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.ticker

    usage = """
    Filters the data based on the threshold for
    similarity and saves a new file which
    contains only the proteases that where under 
    that threshold. 
    """
    print(usage)

    main()
