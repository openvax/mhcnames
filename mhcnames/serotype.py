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

from .gene import Gene


class Serotype(Gene):
    def __init__(self, species_prefix, name, alleles):
        if len(alleles) == 0:
            raise ValueError("Cannot create Serotype without alleles")

        gene_names = {allele.gene_name for allele in alleles}
        if len(gene_names) != 1:
            raise ValueError(
                "Serotype cannot span multiple genes: %s" % (
                    gene_names,))
        gene_name = list(gene_names)[0]
        Gene.__init__(self, species_prefix, gene_name)
        self.name = name
        self.alleles = alleles

    @classmethod
    def field_names(cls):
        return ("species_prefix", "name", "alleles")

    def normalized_string(self, include_species=True):
        if include_species:
            return "%s-%s" % (self.species_prefix, self.name)
        else:
            return self.name

    def compact_string(self, include_species=False):
        return self.normalized_string(include_species=include_species)

    def to_record(self):
        d = Gene.to_record(self)
        d["serotype"] = self.normalized_string()
        d["alleles_in_serotype"] = ";".join([
            allele.normalized_string() for allele in self.alleles
        ])
        return d
