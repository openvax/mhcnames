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
from .allele_name import parse_allele_name, AlleleName
from .allele_parse_error import AlleleParseError

def infer_alpha_chain(beta):
    """
    Given a parsed beta chain of a class II MHC, infer the most frequent
    corresponding alpha chain.
    """
    if beta.gene.startswith("DRB"):
        return AlleleName(species="HLA", gene="DRA1", allele_family="01", allele_code="01")
    elif beta.gene.startswith("DPB"):
        # Most common alpha chain for DP is DPA*01:03 but we really
        # need to change this logic to use a lookup table of pairwise
        # frequencies for inferring the alpha-beta pairing
        return AlleleName(
            species="HLA", gene="DPA1", allele_family="01", allele_code="03")
    elif beta.gene.startswith("DQB"):
        # Most common DQ alpha (according to wikipedia)
        # DQA1*01:02
        return AlleleName(
            species="HLA", gene="DQA1", allele_family="01", allele_code="02")
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