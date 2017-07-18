# Copyright (c) 2017. Mount Sinai School of Medicine
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

from .allele_name import AlleleName
from .class2 import parse_classi_or_classii_allele_name

_normalized_allele_cache = {}

def normalize_allele_name(raw_allele):
    """MHC alleles are named with a frustratingly loose system. It's not uncommon
    to see dozens of different forms for the same allele.

    Note: this function works with both class I and class II allele names (including
    alpha/beta pairs).

    For example, these all refer to the same MHC sequence:
        - HLA-A*02:01
        - HLA-A02:01
        - HLA-A:02:01
        - HLA-A0201
        - HLA-A00201

    Additionally, for human alleles, the species prefix is often omitted:
        - A*02:01
        - A*00201
        - A*0201
        - A02:01
        - A:02:01
        - A:002:01
        - A0201
        - A00201

    We might also encounter "6 digit" and "8 digit" MHC types (which specify
    variants that don't affect amino acid sequence), for our purposes these
    should be truncated to their "4-digit" forms:
        - A*02:01:01
        - A*02:01:01:01
    There are also suffixes which we're going to ignore:
        - HLA-A*02:01:01G

    And lastly, for human alleles, there are serotypes which we'll treat
    as approximately equal to a 4-digit type.
        - HLA-A2
        - A2

    These should all be normalized to:
        HLA-A*02:01
    """
    if raw_allele in _normalized_allele_cache:
        return _normalized_allele_cache[raw_allele]
    parsed_alleles = parse_classi_or_classii_allele_name(raw_allele)
    species = parsed_alleles[0].species
    normalized_list = [species]
    for parsed_allele in parsed_alleles:
        if len(parsed_allele.allele_family) > 0:
            normalized_list.append("%s*%s:%s" % (
                parsed_allele.gene,
                parsed_allele.allele_family,
                parsed_allele.allele_code))
        else:
            # mice don't have allele families
            # e.g. H-2-Kd
            # species = H-2
            # gene = K
            # allele = d
            normalized_list.append("%s%s" % (
                parsed_allele.gene,
                parsed_allele.allele_code))
    normalized = "-".join(normalized_list)
    _normalized_allele_cache[raw_allele] = normalized
    return normalized

_DRA1_0101 = AlleleName(
    species="HLA",
    gene="DRA1",
    allele_family="01",
    allele_code="01")

def compact_allele_name(raw_allele):
    """
    Turn HLA-A*02:01 into A0201 or H-2-D-b into H-2Db or
    HLA-DPA1*01:05-DPB1*100:01 into DPA10105-DPB110001
    """
    parsed_alleles = parse_classi_or_classii_allele_name(raw_allele)
    normalized_list = []
    if len(parsed_alleles) == 2:
        alpha, beta = parsed_alleles
        # by convention the alpha allelle is omitted since it's assumed
        # to be DRA1*01:01
        if alpha == _DRA1_0101:
            parsed_alleles = [beta]

    for parsed_allele in parsed_alleles:
        if len(parsed_allele.allele_family) > 0:
            normalized_list.append("%s%s%s" % (
                parsed_allele.gene,
                parsed_allele.allele_family,
                parsed_allele.allele_code))
        else:
            # mice don't have allele families
            normalized_list.append("%s%s" % (
                parsed_allele.gene,
                parsed_allele.allele_code))
    return "-".join(normalized_list)
