import pandas as pd
from mhcnames import normalized_string, AlleleParseError


def test_parse_all_IEDB_allele_names():
    df = pd.read_csv("iedb_allele_counts.csv")
    n_wrong = 0
    n_total = 0
    n_wrong_iedb_count = 0
    n_total_iedb_count = 0
    for (allele, count) in zip(df.allele, df.number_of_entries):
        n_total += 1
        n_total_iedb_count += count
        try:
            normalized_string(allele)
        except AlleleParseError as e:
            print(e)
            n_wrong += 1
            n_wrong_iedb_count += count
    n_correct = n_total - n_wrong
    n_correct_count = n_total_iedb_count - n_wrong_iedb_count

    print("Parsed %d/%d IEDB allele names, %d/%d (%0.2f%%) of entries" % (
        n_correct,
        n_total,
        n_correct_count,
        n_total_iedb_count,
        100.0 * n_correct_count / n_total_iedb_count))
    assert n_correct == n_total

if __name__ == "__main__":
    test_parse_all_IEDB_allele_names()
