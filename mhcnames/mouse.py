# Copyright (c) 2016. Mount Sinai School of Medicine
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


from .allele_parse_error import AlleleParseError
from .parsing_helpers import parse_separator

def parse_mouse_allele_name(name):
    """Parses mouse MHc alleles such as H2-Kd, H-2-Db, H2-IAb.
    Returns pair of (gene, allele_code).
    """
    original = name
    if name.upper().startswith("H2"):
        name = name[2:]
    elif name.upper().startswith("H-2"):
        name = name[3:]
    _, name = parse_separator(name)

    # special logic for mouse alleles
    if name.upper().startswith("I"):
        # class II mouse allele
        if len(name) < 2:
            raise AlleleParseError("Incomplete mouse MHC allele: %s" % original)
        gene_name = name[:2]
        name = name[2:]
    else:
        # class I mouse allele
        if len(name) < 1:
            raise AlleleParseError("Incomplete mouse MHC allele: %s" % original)
        gene_name = name[0]
        name = name[1:]
    _, name = parse_separator(name)

    if len(name) != 1:
        raise AlleleParseError(
            "Malformed mouse MHC allele: %s, parse error at %s" % (
                original, name))
    allele = name[0]
    return gene_name.upper(), allele.lower()
