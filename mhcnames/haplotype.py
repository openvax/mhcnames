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

from serializable import Serializable

from .mhc_class_restriction import (
    is_valid_restriction,
    restrict_alleles,
    valid_class_types_and_subtypes
)

class Haplotype(Serializable):
    def __init__(self, species_prefix, haplotype_name, alleles, class_restriction=None):
        self.species_prefix = species_prefix
        self.haplotype_name = haplotype_name
        self.alleles = alleles
        assert (
            (class_restriction is None) or
            (class_restriction in valid_class_types_and_subtypes)
        )

    def restrict_alleles(self, class_restriction):
        assert class_restriction is not None
        if self.class_restriction == class_restriction:
            return self
        if not is_valid_restriction(self.class_restriction, class_restriction):
            raise ValueError(
                "Cannot restrict '%s' to class '%s'" % (
                    self.normalized_string(),
                    class_restriction))
        return restrict_alleles(self.alleles, class_restriction)

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
