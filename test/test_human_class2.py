from nose.tools import eq_

from mhcnames import (
    normalized_string,
    compact_string,
)


def test_human_class2():
    expected = "HLA-DRA1*01:01-DRB1*01:02"
    expected_compact = "DRB10102"
    for name in ["DRB1_0102",
                 "DRB101:02",
                 "HLA-DRB1_0102",
                 "DRB10102",
                 "DRB1*0102",
                 "HLA-DRB1*0102",
                 "HLA-DRB1*01:02",
                 "DRB0102"]:
        eq_(normalized_string(name), expected)
        eq_(compact_string(name), expected_compact)


def test_human_class2_alpha_beta():
    expected = "HLA-DPA1*01:05-DPB1*100:01"
    expected_compact = "DPA10105-DPB110001"
    for name in ["DPA10105-DPB110001",
                 "HLA-DPA1*01:05-DPB1*100:01",
                 "hla-dpa1*0105-dpb1*10001",
                 "dpa1*0105-dpb1*10001",
                 "HLA-DPA1*01:05/DPB1*100:01",
                 "DPA10105/DPB110001"]:
        eq_(normalized_string(name), expected)
        eq_(compact_string(name), expected_compact)


def test_alpha_chain_inference_DP():
    # example of common DP haplotype from wikipedia article on HLA-DP
    just_beta = "HLA-DPB1*04:01"
    expected = "HLA-DPA1*01:03-DPB1*04:01"
    eq_(normalized_string(just_beta, infer_class2_pairing=True), expected)
    eq_(normalized_string(just_beta, infer_class2_pairing=False), just_beta)


def test_alpha_chain_inference_DQ():
    # example of common DQ haplotype from wikipedia article on HLA-DQ
    just_beta = "HLA-DQB1*06:02"
    expected = "HLA-DQA1*01:02-DQB1*06:02"
    eq_(normalized_string(just_beta, infer_class2_pairing=True), expected)
    eq_(normalized_string(just_beta, infer_class2_pairing=False), just_beta)
