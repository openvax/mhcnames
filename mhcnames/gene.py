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

from pytypes import typechecked

from .mhc_class_helpers import is_class1, is_class2
from .parsed_result import ParsedResult
from .species import Species


class Gene(ParsedResult):
    @typechecked
    def __init__(self, species : Species, gene_name : str):
        self.species = species
        self.gene_name = gene_name

    @classmethod
    def field_names(cls):
        """
        Returns name of fields used in constructor
        """
        return ("species", "gene_name")

    @property
    def name(self):
        # alias for gene_name
        return self.gene_name

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def mhc_class(self):
        return self.species.get_mhc_class_of_gene(self.gene_name)

    @property
    def is_class1(self):
        return is_class1(self.mhc_class)

    @property
    def is_class2(self):
        return is_class2(self.mhc_class)

    @classmethod
    def get(cls, species_prefix, gene_name : str):
        """
        Returns Gene if gene name is in ontology, otherwise None
        """
        # use the canonical gene name e.g. "A" and not "a"
        if species_prefix.__class__ is Species:
            species = species_prefix
        else:
            species = Species.get(species_prefix)

        if species is None:
            return None
        gene_name = species.find_matching_gene_name(gene_name)
        if gene_name is None:
            return None
        return Gene(species, gene_name)

    def normalized_string(
            self,
            include_species=True,
            use_species_alias=True):
        if include_species:
            # non-classical genes located outside the MHC locus
            # get identified with the common species name if possible
            if self.mhc_class == "Id":
                species_str = self.common_species_name
            elif use_species_alias:
                species_str = self.species.historic_alias
            else:
                species_str = self.species_prefix
            return "%s-%s" % (species_str, self.gene_name)
        else:
            return self.gene_name

    def compact_string(
            self,
            include_species=False,
            use_species_alias=True):
        """
        Compact representation of a Locus, currently same as the
        normalized representation.
        """
        return Gene.normalized_string(
            self,
            include_species=include_species,
            use_species_alias=use_species_alias)


    def to_record(self):
        """
        Returns a user-viewable ordered dictionary with a representation  of
        this gene, and the values of the following methods:
            - is_mutant
            - get_mhc_class
        """
        d = self.species.to_record()
        d["gene"] = self.normalized_string()
        d["mhc_class"] = self.mhc_class
        return d

