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

from .data import load
from .normalizing_dictionary import NormalizingDictionary

# dictionary mapping group -> scientific name -> info dictionary
raw_species_dict = load("species.yaml")

common_names_to_scientific_names = NormalizingDictionary()
prefix_to_scientific_name = NormalizingDictionary()
prefix_to_alias = NormalizingDictionary()
alias_to_four_letter_code = NormalizingDictionary()

for species_group, species_dicts in raw_species_dict.items():
    for scientific_name, species_dict in species_dicts.items():
        prefix = species_dict["prefix"]
        alias = species_dict.get("alias")
        if alias:
            prefixes = [prefix, alias]
            alias_to_four_letter_code[alias] = prefix
            prefix_to_alias[prefix] = alias
        else:
            prefixes = [prefix]

        common_names = species_dict["common name"]
        if type(common_names) is not list:
            common_names = [common_names]
        for common_name in common_names:
            common_names_to_scientific_names[common_name] = scientific_name

        for prefix in prefixes:
            prefix_to_scientific_name[prefix] = scientific_name


scientific_name_to_common_names = common_names_to_scientific_names.invert()

prefix_to_common_names = NormalizingDictionary()
for prefix, scientific_name in prefix_to_scientific_name.items():
    common_names = scientific_name_to_common_names[scientific_name]
    scientific_name_to_common_names[scientific_name] = common_names
    prefix_to_common_names[prefix] = common_names

scientific_name_to_prefixes = prefix_to_scientific_name.invert()
common_name_to_mhc_prefixes = prefix_to_common_names.invert()

# since each species can have more than one common name, it's useful to
# map the MHC prefixes and scientific names to shortest common name


def keep_shortest_name(names):
    return min(names, key=len)

prefix_to_default_common_name = \
    prefix_to_common_names.map_values(keep_shortest_name)
scientific_to_default_common_name = \
    scientific_name_to_common_names.map_values(keep_shortest_name)


# collect all unique scientific, common, and MHC prefix names for species
all_scientific_names = set(scientific_name_to_common_names.keys())
all_prefixes = set(prefix_to_scientific_name.keys())
all_common_names = set(common_names_to_scientific_names.keys())
