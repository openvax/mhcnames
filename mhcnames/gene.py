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

from collections import OrderedDict

from serializable import Serializable

from .species_registry import find_matching_species_info


class Gene(Serializable):
    def __init__(self, species_prefix, gene_name):
        self.species_prefix = species_prefix
        self.gene_name = gene_name

    @property
    def species_info(self):
        return find_matching_species_info(self.species_prefix)

    def get_mhc_class(self):
        return self.species_info.get_mhc_class(self.gene_name)

    def normalized_string(self, include_species=True):
        if include_species:
            return "%s-%s" % (self.species_prefix, self.gene_name)
        else:
            return self.gene_name

    def compact_string(self, include_species=True):
        """
        Compact representation of a Locus, currently same as the
        normalized representation.
        """
        return self.normalized_string(include_species=include_species)

    @classmethod
    def from_tuple(cls, t):
        return cls(
            species_prefix=t[0],
            gene_name=t[1]
        )

    def to_tuple(self):
        return (
            self.species_prefix,
            self.gene_name,
        )

    def to_gene(self):
        """
        Descendant classes use this method to project their fields down
        to a Gene object.
        """
        if self.__class__ is Gene:
            return self
        else:
            return Gene(self.species_prefix, self.gene_name)

    def to_dict(self):
        """
        Returns dictionary with all fields of this gene, its normalized,
        representation and the values of the following methods:
            - is_mutant
            - get_mhc_class
        """
        return OrderedDict([
            ("gene", self.normalized_string()),
            ("mhc_class", self.get_mhc_class()),
            ("is_mutant", False),
        ])
