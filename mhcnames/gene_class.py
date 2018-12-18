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

from .species import Species
from .mhc_class import (
    normalize_mhc_class_string,
    is_class1,
    is_class2,
    is_valid_restriction
)


class GeneClass(Species):
    """
    Wrapper class for species combined with MHC classes such as
    "I", "Ia", "Ib", "II", "IIa", &c
    which provides some utility functions.
    """
    def __init__(self, species_prefix, mhc_class):
        Species.__init__(self, species_prefix)
        self.mhc_class = normalize_mhc_class_string(mhc_class)

    @classmethod
    def field_names(cls):
        return ("species_prefix", "mhc_class")

    @property
    def is_class1(self):
        return is_class1(self.mhc_class)

    @property
    def is_class2(self):
        return is_class2(self.mhc_class)

    def to_dict(self):
        return {
            "species_prefix": self.species_prefix,
            "mhc_class": self.mhc_class,
        }

    def genes(self):
        """
        Returns all gene names whose MHC class matches this MHC class
        """
        return [
            g
            for g in Species.genes(self)
            if is_valid_restriction(self.mhc_class, self.get_mhc_class_of_gene(g))
        ]

    @classmethod
    def from_dict(cls, d):
        return GeneClass(**d)

    def to_record(self):
        d = Species.to_record()
        d["mhc_class"] = self.mhc_class
        return d

    def normalized_string(self, include_species=True):
        if include_species:
            if self.common_species_name:
                species_str = self.common_species_name
            else:
                species_str = self.species_prefix
            return "%s class %s" % (species_str, self.mhc_class)
        else:
            return "class %s" % self.mhc_class

    def compact_string(self, include_species=True):
        """
        Compact representation of an MHC class, currently same as the
        normalized representation.
        """
        return self.normalized_string(include_species=include_species)