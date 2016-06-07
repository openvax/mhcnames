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

from collections import namedtuple

from .parsing_helpers import (
    AlleleParseError,
    parse_separator,
    parse_alphanum,
    parse_numbers,
    parse_letters
)
from .mouse import parse_mouse_allele_name
from .species import split_species_prefix

AlleleName = namedtuple("AlleleName", [
    "species",
    "gene",
    "allele_family",
    "allele_code"
])


def parse_classi_or_classii_allele_name(name):
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
    name = name.replace("_", "")

    parts = name.split("-")
    assert len(parts) <= 2, "Allele has too many parts: %s : %s" % (
        name, parts)
    if len(parts) == 1:
        return (parse_allele_name(name, species),)
    else:
        return (parse_allele_name(parts[0], species),
                parse_allele_name(parts[1], species))


def parse_allele_name(name, species_prefix=None):
    """Takes an allele name and splits it into four parts:
        1) species prefix
        2) gene name
        3) allele family
        4) allele code

    If species_prefix is provided, that is used instead of getting the species prefix from the name.
    (And in that case, a species prefix in the name will result in an error being raised.)

    For example, in all of the following inputs:
        "HLA-A*02:01"
        "A0201"
        "A00201"
    The result is a AlleleName object. Example:
        AlleleName(
            species="HLA",  # species prefix
            gene="A",    # gene name
            allele_family="02",   # allele family
            allele_code="01",   # allele code
        )

    The logic for other species mostly resembles the naming system for humans,
    except for mice, rats, and swine, which have archaic nomenclature.
    """
    original = name
    name = name.strip()

    if len(name) == 0:
        raise ValueError("Can't normalize empty MHC allele name")

    species_from_name, name = split_species_prefix(name)
    if species_prefix:
        if species_from_name:
            raise ValueError("If a species is passed in, we better not have another "
                             "species in the name itself.")
        species = species_prefix
    else:
        species = species_from_name

    if species in ("H-2", "H2"):
        gene, allele_code = parse_mouse_allele_name("H-2-" + name)
        # mice don't have allele families
        return AlleleName("H-2", gene, "", allele_code)

    if len(name) == 0:
        raise AlleleParseError("Incomplete MHC allele name: %s" % (original,))

    elif not species:
        # assume that a missing species name means we're dealing with a
        # human HLA allele
        if "-" in name:
            raise AlleleParseError("Can't parse allele name: %s" % original)
        species = "HLA"

    if name[0].upper() == "D":
        # MHC class II genes like "DQA1" need to be parsed with both
        # letters and numbers
        gene, name = parse_alphanum(name, 4)
    elif name[0].isalpha():
        # if there are more separators to come, then assume the gene names
        # can have the form "DQA1"
        gene, name = parse_letters(name)
    elif name[0].isdigit():
        gene, name = parse_numbers(name)
    else:
        raise AlleleParseError(
            "Can't parse gene name from allele: %s" % original)

    if len(gene) == 0:
        raise AlleleParseError("No MHC gene name given in %s" % original)
    if len(name) == 0:
        raise AlleleParseError("Malformed MHC type %s" % original)

    gene = gene.upper()

    # skip initial separator
    sep, name = parse_separator(name)
    # if all that's left is e.g. "0201" then only parse the
    # first two digits as the family code
    if len(name) == 4:
        family, name = parse_numbers(name, max_len=2)
    else:
        family, name = parse_numbers(name, max_len=3)

    sep, name = parse_separator(name)

    allele_code, name = parse_numbers(name)

    if len(family) == 1:
        family = "0" + family
    elif len(family) == 3 and family[0] == "0":
        family = family[1:]

    if len(allele_code) == 0:
        allele_code = "01"
    elif len(allele_code) == 1:
        # change HLA-A*2:01 into HLA-A*02:01
        allele_code = "0" + allele_code
    elif len(allele_code) == 3 and allele_code[0] == "0":
        # normalize HLA-A*002:01 into HLA-A*02:01
        allele_code = allele_code[1:]
    return AlleleName(species, gene, family, allele_code)

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

def compact_allele_name(raw_allele):
    """
    Turn HLA-A*02:01 into A0201 or H-2-D-b into H-2Db or
    HLA-DPA1*01:05-DPB1*100:01 into DPA10105-DPB110001
    """
    parsed_alleles = parse_classi_or_classii_allele_name(raw_allele)
    normalized_list = []
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
