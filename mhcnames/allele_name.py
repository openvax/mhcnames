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
    parse_separator,
    parse_alphanum,
    parse_numbers,
    parse_letters
)
from .mouse import parse_mouse_allele_name
from .species import split_species_prefix
from .allele_parse_error import AlleleParseError

AlleleName = namedtuple("AlleleName", [
    "species",
    "gene",
    "allele_family",
    "allele_code"
])


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
        # TODO: make a list of known species/gene pairs, along with
        # gene synonyms. That should significantly imporve on this kind of
        # ad-hoc synonym handling.
        if gene == "DRA":
            gene = "DRA1"
        elif gene == "DRB":
            gene = "DRB1"
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
    # If all that's left is e.g. "0201" then only parse the
    # first two digits as the family code. Also, human Class I alleles
    # seem to be exceptional in that they have only 2 digit allele
    # families but 3 digit allele codes
    # (other species like sheep have 3 digits followed by 2 digits)
    if len(name) == 4 or (species == "HLA" and gene in ("A", "B", "C")):
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
        # normalize HLA-A*02:001 into HLA-A*02:01
        allele_code = allele_code[1:]

    return AlleleName(species, gene, family, allele_code)
