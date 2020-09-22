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

from .normalizing_dictionary import NormalizingDictionary
from .data import serotypes


def _create_serotype_aliases_dict():
    species_to_alias_to_serotype = NormalizingDictionary()
    for (species, species_serotypes) in serotypes.items():
        species_to_alias_to_serotype[species] = NormalizingDictionary()
        for serotype_name in species_serotypes:
            aliases = {serotype_name}
            if serotype_name.startswith("Cw"):
                aliases.add("C" + serotype_name[2:])
            if serotype_name.startswith("DPw"):
                aliases.add("DP" + serotype_name[3:])
            aliases.add(serotype_name.upper())
            for alias in aliases:
                species_to_alias_to_serotype[species][alias] = serotype_name
    return species_to_alias_to_serotype

species_prefix_to_alias_to_serotype = _create_serotype_aliases_dict()


def get_serotype(species_prefix, serotype_name):
    """
    Returns either None (if serotype doesn't exist) or tuple with following
    entries:
        - normalized species prefix
        - normalized serotype name
        - list of alleles in serotype
    """
    alias_to_serotype_dict = species_prefix_to_alias_to_serotype.get(species_prefix, {})
    corrected_serotype = alias_to_serotype_dict.get(serotype_name)
    if corrected_serotype is None:
        return None
    else:
        corrected_species_prefix = \
            species_prefix_to_alias_to_serotype.original_key(species_prefix)
        allele_list = serotypes[corrected_species_prefix][corrected_serotype]
        return (
            corrected_species_prefix,
            corrected_serotype,
            allele_list
        )
