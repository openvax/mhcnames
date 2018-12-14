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
from .parsed_result import ParsedResult


class AlphaBetaPair(ParsedResult):
    def __init__(self, alpha, beta):
        self.alpha = alpha
        self.beta = beta

    @classmethod
    def field_names(cls):
        return ("alpha", "beta")

    @property
    def species_prefix(self):
        return self.alpha.species_prefix

    def normalized_string(self, include_species=True):
        return "%s/%s" % (
            self.alpha.normalized_string(include_species=include_species),
            self.beta.normalized_string(include_species=False))

    def compact_string(self, include_species=False):
        return "%s/%s" % (
            self.alpha.compact_string(include_species=include_species),
            self.beta.compact_string(include_species=False))

    @property
    def mhc_class(self):
        alpha_class = self.alpha.mhc_class
        if alpha_class:
            return alpha_class
        else:
            return self.beta.mhc_class

    def to_record(self):
        # return a dictionary that has the same elements as Gene.to_dict()
        # along with "allele"
        return OrderedDict([
            ("gene", "%s/%s" % (self.alpha.gene_name, self.beta.gene_name)),
            ("mhc_class", self.mhc_class),
            ("is_mutant", False),
            ("allele", self.normalized_string()),
        ])
