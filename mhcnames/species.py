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

from .data import gene_ontology as gene_ontology_dict

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


def create_aliases():
    aliases = exemplar_species.copy()

    for key, value in list(gene_ontology_dict.items()):
        upper = key.upper()
        if upper != key:
            aliases[upper] = key
        lower = key.lower()
        if lower != key:
            aliases[lower] = key
    return aliases

species_aliases = create_aliases()

def find_matching_species_prefix(name):
    """
    Returns either normalized species prefix or None
    """
    if name in gene_ontology_dict:
        return name

    upper_no_dash = name.upper().replace("-", "")

    if upper_no_dash in species_aliases:
        return gene_name_aliases[upper_no_dash]

    return None

def get_species_info(name):
    """
    Returns a dictionary with gene lists for classes "Ia", "Ib", "IIa", "IIb"
    or None if species can't be found
    """
    # change H-2 -> H2 and RT-1 to RT1
    species_prefix = find_matching_species_prefix(name)
    if species_prefix:
        return gene_ontology_dict[species_prefix]
    else:
        raise ValueError("Could not find species information for '%s'" % name)

def get_mhc_class(species_name, gene_name):
    """
    Returns either one of "Ia", "Ib", "IIa", "IIb" or None
    if species can't be found
    """
    species_info = get_species_info(species_name)
    for mhc_class, mhc_class_members in species_info.items():
        if mhc_class == "Ia" or mhc_class == "Ib":
            for member_gene in mhc_class_members:
                if member_gene == gene_name:
                    return mhc_class
        elif mhc_class == "IIa" or mhc_class == "IIb":
            for locus, genes in mhc_class_members.items():
                if locus == gene_name:
                    return mhc_class
                for candidate_gene_name in genes:
                    if candidate_gene_name == gene_name:
                        return mhc_class
    return None

def infer_species_prefix(allele_string):
    """
    Example the beginning of string to see if it matches any of the
    known species prefixes.

    Examples:
        * "hla-a" -> "HLA"
        * "RT1" -> "RT1"
        * "RT-1" -> "RT1"
        * "rano" -> "Rano"
        * "H-2" -> "H2"

    """
    for n in [2, 3, 4]:
        substring = allele_string[:n]
        upper_no_dash = substring.upper().replace("-", "")
        if upper_no_dash in special_prefixes:
            return upper_no_dash
        species_prefix = find_matching_species_prefix(upper_no_dash)
        if species_prefix is not None:
            return species_prefix
    return None
