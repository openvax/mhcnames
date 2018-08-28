# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

from .species import split_species_prefix
from .allele_parse_error import AlleleParseError
from .allele_modifiers import valid_allele_modifers

import re

default_alpha_chains = {
    # The DR alpha chain is effectively constant across the population
    "DRB": "HLA-DRA1*01:01",
    # Most common alpha chain for DP is DPA*01:03 but we really
    # need to change this logic to use a lookup table of pairwise
    # frequencies for inferring the alpha-beta pairing
    "DPB": "HLA-DPA1*01:03",
    # Most common DQ alpha (according to wikipedia) is DQA1*01:02
    # but like DPA we should use pair frequencies in the future
    "DQB": "HLA-DQA1*01:02",
}


def infer_alpha_chain(beta):
    """
    Given a parsed beta chain of a class II MHC, infer the most frequent
    corresponding alpha chain.
    """
    for (gene, alpha) in default_alpha_chains.items():
        if beta.startswith(gene):
            return alpha
    return None

def parse_classi_or_classii_allele_name(name, infer_pair=True):
    """
    Handle different forms of both single and alpha-beta allele names.
    Alpha-beta alleles may look like:

    DPA10105-DPB110001
    HLA-DPA1*01:05-DPB1*100:01
    hla-dpa1*0105-dpb1*10001
    dpa1*0105-dpb1*10001
    HLA-DPA1*01:05/DPB1*100:01

    Other class II alleles may look like:

    DRB1_0102
    DRB101:02
    HLA-DRB1_0102
    """
    name = name.upper().strip()

    # strip off leading prefix info
    prefix_match = re.match("^HLA[-\:\*_]+", name)
    if prefix_match is not None:
        _, prefix_end_offset = prefix_match.span()
        name = name[prefix_end_offset:]

    # gene should start with a letter and can optionally end in a number
    gene_regex = "([A-Za-z][A-Za-z]*[0-9]?)"
    # alleles are expected to be numbers
    allele_regex = "([0-9]+)"
    modifier_regex = "([" + valid_allele_modifers + "])?"
    sep_regex = "[-_\:\*]?"

    two_digit_regex = gene_regex + sep_regex + allele_regex
    four_digit_regex = (
        gene_regex + sep_regex +
        allele_regex + sep_regex +
        allele_regex + modifier_regex
    )
    six_digit_regex = (
        gene_regex + sep_regex +
        allele_regex + sep_regex +
        allele_regex + sep_regex +
        allele_regex + sep_regex + modifier_regex)
    eight_digit_regex = (
        gene_regex + sep_regex +
        allele_regex + sep_regex +
        allele_regex + sep_regex +
        allele_regex + sep_regex +
        allele_regex + sep_regex + modifier_regex
    )
    two_digit_paired_regex = two_digit_regex + sep_regex + two_digit_regex
    four_digit_paired_regex = four_digit_regex + sep_regex + four_digit_regex
    six_digit_paired_regex = six_digit_regex + sep_regex + six_digit_regex
    eight_digit_paired_regex = eight_digit_regex + sep_regex + eight_digit_regex

    eight_digit_paired_match = re.full_match(eight_digit_paired_regex, name):
    if eight_digit_paired_match:
        alpha_gene = eight_digit_paired_match.groups(1)
        alpha_locus = eight_digit_paired_match.groups(2)
        alpha_protein = eight_digit_paired_match.groups(3)
        alpha_coding_sequence = eight_digit_paired_match.groups(4)
        alpha_genomic_sequence = eight_digit_paired_match.groups(5)

        beta_gene = eight_digit_paired_match.groups(6)
        beta_locus = eight_digit_paired_match.groups(7)
        beta_protein = eight_digit_paired_match.groups(8)
        beta_coding_sequence = eight_digit_paired_match.groups(9)
        beta_genomic_sequence = eight_digit_paired_match.groups(10)

    species, name = split_species_prefix(name)

    # Handle the case where alpha/beta pairs are separated with a /.
    name = name.replace("/", "-")

    # Ignored underscores, such as with DRB1_0102
    name = name.replace("_", "*")

    parts = name.split("-")

    if len(parts) == 2:
        alpha_string, beta_string = parts
        alpha = parse_allele_name(alpha_string)
        beta = parse_allele_name(beta_string)
        return (alpha, beta)
    elif len(parts) == 1:
        parsed = parse_allele_name(name, species)
        if parsed.species == "HLA" and infer_pair:
            alpha = infer_alpha_chain(parsed)
            if alpha is not None:
                return (alpha, parsed)
        return (parsed,)
    else:
        raise AlleleParseError(
            "Allele has too many parts: %s" % name)
