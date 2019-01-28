# Copyright (c) 2018-2019. Mount Sinai School of Medicine
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

from .gene import Gene


class NamedAllele(Gene):
    """
    Some species, such as mouse (H2) and rats (RT1) do not yet use
    the standard nomenclature format (e.g. Species-Gene*Group:Protein).

    Instead they identify each distinct MHC molecule with a letter, or sometimes
    short sequence. Example: H2-Kk, RT1-9.5*f

    Also, some older swine (SLA) alleles seem to have not been translated into
    the updated nomenclature (e.g. SLA-1-CHANGDA)
    """
    def __init__(self, species_prefix, gene_name, allele_name):
        Gene.__init__(self, species_prefix, gene_name)
        self.allele_name = allele_name

    def field_names(self):
        return ("species_prefix", "gene_name", "allele_name")

    def to_tuple(self):
        return (self.species_prefix, self.gene_name, self.allele_name)

    @classmethod
    def from_tuple(cls, t):
        return cls(t[0], t[1], t[2])

    def normalized_string(self, include_species=True):
        return "%s%s" % (
            Gene.normalized_string(
                self,
                include_species=include_species),
            self.allele_name)

    def compact_string(self, include_species=False):
        return "%s%s" % (
            Gene.compact_string(
                self,
                include_species=include_species),
            self.allele_name)

    def to_record(self):
        d = Gene.to_record(self)
        d["allele"] = self.normalized_string()
