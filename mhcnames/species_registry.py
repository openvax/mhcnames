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

from .data import gene_ontology as raw_gene_ontology_dict

from .species import Species

# Many old-fasioned naming systems like "equine" ELA now correspond
# to multiple species. For each species-ambiguous prefix, map it to the
# species which has the most complete gene annotations.
exemplar_species = {
    "DLA": "Calu",
    "ELA": "Eqca",
    "OLA": "Ovar",
    "SLA": "Susc",
    "RT1": "Rano"
}

special_prefixes = set(exemplar_species.keys())

# map prefixes to Species objects
species_dict = {}
for species_prefix, species_gene_ontology in raw_gene_ontology_dict.items():
    species_dict[species_prefix] = Species(species_prefix)


def create_species_aliases():
    """
    Create dictionary of species aliases from both the species->MHC class->gene
    ontology and the set of MHC names (like RT1) which actually represent
    multiple species and require an exemplar to be chosen for its gene
    metadata.
    """
    aliases = species_dict.copy()
    for group_name, species_name in exemplar_species.items():
        aliases[group_name] = species_dict[species_name]

    for key, species in list(species_dict.items()):
        upper = key.upper()
        if upper != key:
            aliases[upper] = species
        upper_no_dash = upper.replace("-", "")
        if upper_no_dash not in {key, upper}:
            aliases[upper_no_dash] = species
    return aliases

# dictionary mapping alias names to Species objects
species_aliases_dict = create_species_aliases()


def find_matching_species(name):
    """
    Returns either SpeciesInfo object or None if species can't be found
    """
    if name in species_dict:
        return species_dict[name]
    return species_aliases_dict.get(name.upper().replace("-", ""))


def find_matching_species_prefix(name):
    """
    Returns normalized prefix for given species but will keep
    RT1, DLA, ELA, OLA, SLA from normalizing to a particular member
    of their species group.
    """
    upper_no_dash = name.upper().replace("-", "")
    if upper_no_dash in exemplar_species:
        return upper_no_dash
    species = find_matching_species(name)

    if species is not None:
        return species.species_prefix
    else:
        return None


def infer_species_prefix_substring(name):
    """
    Trying to parse prefixes of alleles such as:
        HLA-A
    but also ones with dashes in the species prefix:
        H-2-K
    and also those lacking any dashes such as:
        H2K

     ...we also need to consider that alleles, haplotypes, etc may come
     immediately after the gene:
        H2Kk
        HLA-A0201

    Returns the normalized species prefix and the original string that matched
    it or None.
    """
    # Try parsing a few different substrings to get the species,
    # and then use the species gene list to determine what the gene is in this string
    candidate_species_substrings = [name]

    if "-" in name:
        # if name is "H-2-K" then try parsing "H" and "H-2" as a species
        # prefix
        parts_split_by_dash = name.split("-")
        candidate_species_substrings.extend([
            parts_split_by_dash[0],
            parts_split_by_dash[0] + "-" + parts_split_by_dash[1]
        ])
    for seq in candidate_species_substrings:
        for n in [2, 3, 4]:
            original_prefix = seq[:n]
            normalized_prefix = find_matching_species_prefix(name[:n])
            if normalized_prefix is not None:
                return normalized_prefix, original_prefix
    return None
