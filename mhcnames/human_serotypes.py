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

from .four_digit_allele import FourDigitAllele
from .data import human_serotypes as raw_serotypes_dict
from .serotype import Serotype

def create_human_serotypes_dict():
    human_serotypes_dict = {}
    for serotype_name, allele_strings in raw_serotypes_dict.items():
        assert serotype_name not in human_serotypes_dict
        alleles = [
            FourDigitAllele.parse(allele) for allele in allele_strings
        ]
        serotype = Serotype(
            name=serotype_name,
            alleles=alleles)
        human_serotypes_dict[serotype_name] = serotype
    return human_serotypes_dict

human_serotypes_dict = create_human_serotypes_dict()

def get_human_serotype_if_exists(name):
    """
    Try to match the given serotype name with one of the existing names,
    return the Serotype object if it exists.

    Normalize "C7" into Cw7" and "DP1" into "DPw1"
    """
    while name.startswith("HLA-"):
        name = name[4:]
    name = name.upper()
    if name.startswith("CW"):
        name = "Cw" + name[2:]
    elif name.startswith("DPW"):
        name = "DPw" + name[3:]
    elif name[0] == "C" and name[1:].isnumeric():
        name = "Cw" + name[1:]
    elif name[:2] == "DP" and name[2:].isnumeric():
        name = "DPw" + name[2:]
    return human_serotypes_dict.get(name)
