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

from .parsed_result import ParsedResult

from .mhc_class_helpers import (
    is_valid_restriction,
    restrict_alleles,
)


class Haplotype(ParsedResult):
    def __init__(
            self,
            species,
            haplotype_name,
            alleles,
            class_restriction=None):
        self.species = species
        self.haplotype_name = haplotype_name
        self.alleles = alleles
        self.class_restriction = class_restriction

    @classmethod
    def field_names(cls):
        return (
            "species",
            "haplotype_name",
            "alleles",
            "class_restriction")

    @property
    def species_prefix(self):
        return self.species.prefix

    def restrict_mhc_class(self, class_restriction):
        assert class_restriction is not None
        if self.class_restriction == class_restriction:
            return self
        if not is_valid_restriction(self.class_restriction, class_restriction):
            raise ValueError(
                "Cannot restrict '%s' to class '%s'" % (
                    self.normalized_string(),
                    class_restriction))
        restricted_alleles = restrict_alleles(self.alleles, class_restriction)
        return Haplotype(
            self.species_prefix,
            self.haplotype_name,
            restricted_alleles,
            class_restriction)

    def normalized_string(self, include_species=True):
        if include_species:
            result = "%s-%s" % (self.species_prefix, self.haplotype_name)
        else:
            result = self.haplotype_name

        if self.class_restriction:
            result += " class %s" % (self.class_restriction,)

        return result

    def compact_string(self, include_species=False):
        return self.normalized_string(include_species=include_species)
