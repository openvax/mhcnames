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

from .mhc_class_helpers import is_class1, is_class2
from .species import Species


class Gene(Species):
    def __init__(self, species_prefix, gene_name):
        Species.__init__(self, species_prefix)
        self.gene_name = gene_name

    @property
    def mhc_class(self):
        return self.get_mhc_class_of_gene(self.gene_name)

    @property
    def is_class1(self):
        return is_class1(self.mhc_class)

    @property
    def is_class2(self):
        return is_class2(self.mhc_class)

    def normalized_string(self, include_species=True):
        if include_species:
            # non-classical genes located outside the MHC locus
            # get identified with the common species name if possible
            if self.mhc_class == "Id":
                species_str = self.common_species_name
            else:
                species_str = self.species_prefix
            return "%s-%s" % (species_str, self.gene_name)
        else:
            return self.gene_name

    def compact_string(self, include_species=False):
        """
        Compact representation of a Locus, currently same as the
        normalized representation.
        """
        return Gene.normalized_string(self, include_species=include_species)

    @classmethod
    def field_names(cls):
        """
        Returns name of fields used in constructor
        """
        return ("species_prefix", "gene_name")

    def to_gene(self):
        """
        Descendant classes use this method to project their fields down
        to a Gene object.
        """
        if self.__class__ is Gene:
            return self
        else:
            return Gene(self.species_prefix, self.gene_name)

    def to_record(self):
        """
        Returns a user-viewable ordered dictionary with a representation  of
        this gene, and the values of the following methods:
            - is_mutant
            - get_mhc_class
        """
        d = Species.to_record()
        d["gene"] = self.normalized_string()
        d["mhc_class"] = self.mhc_class
        return d
