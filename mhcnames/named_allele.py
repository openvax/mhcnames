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

from .allele_group import AlleleGroup

class NamedAllele(AlleleGroup):

    def __init__(self, species_prefix, gene_name, allele_name):
        AlleleGroup.__init__(self, species_prefix, gene_name)
        self.allele_name = allele_name

    def to_dict(self):
        pass

    def parse_substring(self):
        pass

    def parse(self):
        pass

    def normalized_string(self):
        pass

    def compact_string(self):
        pass
