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

"""
Implements parsing logic for the compact representation of
standard naming.

Difference between compact and standard representation:

Result type           |     Standard        |    Compact
-----------------------------------------------------------
    AlleleGroup       |  HLA-A*02           |  HLA-A02
    FourDigitAllele   |  HLA-A*02:01        |  HLA-A0201
    SixDigitAllele    |  HLA-A*02:01:01     |  HLA-A020101
    EightDigitAllele  |  HLA-A*02:01:01:01  |  HLA-A02010101

The biggest source of ambiguity without ':' digit groups is when
one of them actually have three digits, e.g.
    HLA-DRB1*03033 (parsed as HLA-DRB1*03:33)
    Eqca-1*00101 (parsed as Eqca-1*01:01)
    BoLA-2*01201 (parsed as BoLA-2*12:01)
    DLA-88*50801 (parsed as, I think, DLA-88*508:01) TODO: verify parsed form

To handle this correctly we need a species/gene specific piece of information
regarding which fields (the allele group vs. protein ID) can be 3 digits.

"""

from __future__ import print_function, division, absolute_import
