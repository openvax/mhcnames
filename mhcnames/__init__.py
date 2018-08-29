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

from .normalization import normalized_string, compact_string
from .allele_parse_error import AlleleParseError
from .parsing import parse
from .locus import Locus
from .allele_group import AlleleGroup
from .four_digit_allele import FourDigitAllele
from .six_digit_allele import SixDigitAllele
from .eight_digit_allele import EightDigitAllele
from .named_allele import NamedAllele
from .mutation import Mutation
from .mutant_allele import MutantAllele
from .alpha_beta_pair import AlphaBetaPair

__version__ = "1.0.0"

__all__ = [
    "normalized_string",
    "compact_string",
    "AlleleParseError",
    "parse",
    "Locus",
    "AlleleGroup",
    "FourDigitAllele",
    "SixDigitAllele",
    "EightDigitAllele",
    "AlphaBetaPair",
    "Mutation",
    "MutantAllele",
    "NamedAllele",
]
