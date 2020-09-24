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

from .compat import normalized_string, compact_string
from .parse_error import ParseError
from .parser import Parser, parse
from .gene import Gene
from .numeric_alleles import (
    AlleleGroup,
    TwoDigitAllele,
    FourDigitAllele,
    SixDigitAllele,
    EightDigitAllele
)
from .named_allele import NamedAllele
from .mutation import Mutation
from .mutant_allele import MutantAllele
from .alpha_beta_pair import AlphaBetaPair
from .mhc_class import MhcClass
from .serotype import Serotype
from .species import Species
from .haplotype import Haplotype
from .dataframe import dataframe_from_list

__version__ = "1.0.0"

__all__ = [
    "normalized_string",
    "compact_string",
    "ParseError",
    "parse",
    "Gene",
    "MhcClass",
    "AlleleGroup",
    "TwoDigitAllele",
    "FourDigitAllele",
    "SixDigitAllele",
    "EightDigitAllele",
    "AlphaBetaPair",
    "Mutation",
    "MutantAllele",
    "NamedAllele",
    "Haplotype",
    "Serotype",
    "Species",
    "dataframe_from_list",
    "Parser"
]
