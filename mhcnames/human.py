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

from .four_digit_allele import FourDigitAllele
from .data import serotypes as multispecies_serotypes_dict
from .serotype import Serotype

def create_human_serotypes_dict():
    human_serotypes_dict = {}
    for serotype_name, allele_strings in multispecies_serotypes_dict["HLA"].items():
        assert serotype_name not in human_serotypes_dict
        alleles = [
            FourDigitAllele.parse(allele) for allele in allele_strings
        ]
        serotype = Serotype(
            name=serotype_name,
            alleles=alleles)
        human_serotypes_dict[serotype_name] = serotype
    return human_serotypes_dict

def create_serotype_aliases(serotypes):
    aliases = {}
    actions = [
        # include HLA-CW7 as an alias for the serotype Cw7
        lambda s: s.upper(),
        # include HLA-Cw7 as an alias for the serotype Cw7
        lambda s: "HLA-" + s,
        # include C7 an alias for serotype Cw7
        lambda s: s.replace("W", "").replace("w", ""),
    ]
    original_serotype_set = set(serotypes)
    for fn in actions:
        for original_serotype in original_serotype_set:
            modified = fn(original_serotype)
            if modified not in aliases:
                aliases[modified] = original_serotype
        for (other_alias, original_serotype) in aliases.items():
            modified = fn(other_alias)
            if modified not in aliases:
                aliases[modified] = original_serotype
    return aliases

# contains mapping from serotype name to Serotype object
human_serotypes_dict = create_human_serotypes_dict()
human_serotype_aliases_dict = create_serotype_aliases(human_serotypes_dict)

def get_human_serotype_if_exists(name):
    """
    Try to match the given serotype name with one of the existing names,
    return the Serotype object if it exists.

    Normalize "C7" into Cw7" and "DP1" into "DPw1"
    """
    name = name.upper()
    if name in human_serotype_aliases_dict:
        name = human_serotype_aliases_dict[name]
    return human_serotypes_dict.get(name)

def is_human(name):
    upper = name.upper()
    return upper.startswith("HLA") or get_human_serotype_if_exists(name)
