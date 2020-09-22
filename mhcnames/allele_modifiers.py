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

# N = null
# L = low surface expression
# S = soluble / secreted
# C = cytoplasm / not on the surface
# A = aberrant expression
# Q = questionable
# G = group of alleles with identical peptide binding region
#
# For more information see:
# - http://hla.alleles.org/nomenclature/naming.html
# - http://hla.alleles.org/alleles/g_groups.html

valid_allele_modifiers = "NLSCAQG"

allele_modifier_regex_string = "[%s]" % valid_allele_modifiers
allele_modifier_regex_group_string = "(%s)" % allele_modifier_regex_string


def check_for_allele_modifier(seq):
    """
    See if single remaining character after parsing an allele is one
    of the valid modifier characters.
    """
    if len(seq) == 1:
        uppercase = seq.upper()
        for candidate_modifier in valid_allele_modifiers:
            if candidate_modifier == uppercase:
                return candidate_modifier, ""
    else:
        return None, seq
