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

from __future__ import print_function, division, absolute_import

from collections import defaultdict

from .data import species
from .helpers import (
    expand_with_uppercase_and_no_dash,
    apply_string_expansion_to_set_members
)

species_scientific_names_to_common_names = defaultdict(set)
species_common_names_to_scientific_names = {}
species_scientific_names_to_mhc_prefixes = {}
species_common_names_to_mhc_prefixes = {}
species_mhc_prefixes_to_scientific_names = {}
species_mhc_prefixes_to_common_names = {}
species_mhc_prefixes_to_aliases = {}
species_mhc_prefix_aliases_to_four_letter_codes = {}
all_species_common_names = set([])
all_species_common_names_with_case_variants = set([])

all_species_scientific_names = set([])
all_species_scientific_names_with_case_variants = set([])

all_mhc_prefixes = set([])
all_mhc_prefixes_with_case_variants = set([])

for species_group, species_dicts in species.items():
    for scientific_name, species_dict in species_dicts.items():
        all_species_scientific_names.add(scientific_name)
        prefix = species_dict["prefix"]
        alias = species_dict.get("alias")
        if alias:
            default_prefix = alias
            prefixes = [prefix, alias]
            species_mhc_prefix_aliases_to_four_letter_codes[alias] = prefix
            species_mhc_prefixes_to_aliases[prefix] = alias
        else:
            default_prefix = prefix
            prefixes = [prefix]

        common_names = species_dict["common name"]
        if type(common_names) is not list:
            common_names = [common_names]

        for common_name in common_names:
            all_species_common_names.add(common_names)
            for x in expand_with_uppercase_and_no_dash(common_name):
                species_common_names_to_mhc_prefixes[x] = default_prefix
        default_common_name = common_names[0]
        species_scientific_names_to_common_names[scientific_name] = default_common_name

        for prefix in prefixes:
            all_mhc_prefixes.add(prefix)
            species_scientific_names_to_mhc_prefixes[scientific_name].extend(prefixes)
            species_mhc_prefixes_to_scientific_names[prefix] = scientific_name
            species_mhc_prefixes_to_common_names[prefix] = default_common_name

"""
def find_matching_species_prefix(s):
    if s in
"""