# Copyright (c) 2018. Mount Sinai School of Medicine
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

"""
Implements parsing logic for the newer standard naming for
allele groups (e.g. HLA-A*02), four digit alleles (e.g. HLA-A*02:01),
six digit alleles (e.g. HLA-A*02:01:01) and eight digit alleles
(e.g. HLA-A*02:01:01:01)
"""

from __future__ import print_function, division, absolute_import

import re

from .gene import Gene
from .allele_group import AlleleGroup
from .four_digit_allele import FourDigitAllele
from .six_digit_allele import SixDigitAllele
from .eight_digit_allele import EightDigitAllele
from .allele_modifiers import allele_modifier_regex_group_string


species_prefix_regex_string = "([A-Z][A-Za-z0-9]+)"

gene_regex_string = "([A-Za-z0-9]+\.?-?[0-9]*)"

species_with_gene_regex_string = "%s-%s" % (
    species_prefix_regex_string,
    gene_regex_string)

gene_regex = re.compile(species_with_gene_regex_string)

allele_group_regex_string = species_with_gene_regex_string + "\*(\d+)"

allele_group_regex = re.compile(allele_group_regex_string)

# four digit, six digit, and eight digit alleles can optionally have a modifier
# at the end of the allele spec
optional_modifier = allele_modifier_regex_group_string + "?"

four_digit_regex_string = allele_group_regex_string + ":(\d\d+)"

four_digit_regex_string_with_modifier = four_digit_regex_string + optional_modifier
four_digit_regex = re.compile(four_digit_regex_string_with_modifier)

six_digit_regex_string = four_digit_regex_string + ":(\d\d+)"
six_digit_regex_string_with_modifier = six_digit_regex_string + optional_modifier
six_digit_regex = re.compile(six_digit_regex_string_with_modifier)

eight_digit_regex_string = six_digit_regex_string + ":(\d\d+)"
eight_digit_regex_string_with_modifier = eight_digit_regex_string + optional_modifier
eight_digit_regex = re.compile(eight_digit_regex_string_with_modifier)


def parse_standard_allele_name(name):
    if "-" not in name:
        # if sequence isn't even like "HLA-A"
        return None
    num_star = name.count("*")
    num_colon = name.count(":")

    if num_star == 0 and num_colon == 0:
        regex, result_class = (gene_regex, Gene)
    elif num_star == 1 and num_colon == 0:
        regex, result_class = (allele_group_regex, AlleleGroup)
    elif num_star == 1 and num_colon == 1:
        regex, result_class = (four_digit_regex, FourDigitAllele)
    elif num_star == 1 and num_colon == 2:
        regex, result_class = (six_digit_regex, SixDigitAllele)
    elif num_star == 1 and num_colon == 3:
        regex, result_class = (eight_digit_regex, EightDigitAllele)
    else:
        return None
    match = regex.fullmatch(name)
    print(regex, result_class.__name__, name, match)
    if match:
        return result_class.from_tuple(match.groups())
    return None
