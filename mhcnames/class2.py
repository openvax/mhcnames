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

from .species import split_species_prefix
from .allele_name import parse_allele_name, AlleleName
from .allele_parse_error import AlleleParseError

def parse_classi_or_classii_allele_name(name):
    """
    Handle different forms of both single and alpha-beta allele names.
    Alpha-beta alleles may look like:

    DPA10105-DPB110001
    HLA-DPA1*01:05-DPB1*100:01
    hla-dpa1*0105-dpb1*10001
    dpa1*0105-dpb1*10001
    HLA-DPA1*01:05/DPB1*100:01

    Other class II alleles may look like:

    DRB1_0102
    DRB101:02
    HLA-DRB1_0102
    """
    species, name = split_species_prefix(name)

    # Handle the case where alpha/beta pairs are separated with a /.
    name = name.replace("/", "-")

    # Ignored underscores, such as with DRB1_0102
    name = name.replace("_", "*")

    parts = name.split("-")
    if len(parts) > 2:
        raise AlleleParseError(
            "Allele has too many parts: %s" % name)
    if len(parts) == 1:
        parsed = parse_allele_name(name, species)
        if parsed.species == "HLA" and parsed.gene.startswith("DRB"):
            alpha = AlleleName(
                species="HLA",
                gene="DRA1",
                allele_family="01",
                allele_code="01")
            return (alpha, parsed)
        else:
            return (parsed,)
    else:
        alpha_string, beta_string = parts
        alpha = parse_allele_name(alpha_string)
        beta = parse_allele_name(beta_string)
        return (alpha, beta)
