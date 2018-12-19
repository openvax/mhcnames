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

import pandas as pd

from .parsing import parse


def normalized_string(
        raw_string,
        include_species_prefix=True,
        infer_class2_pairing=True,
        default_species_prefix="HLA"):
    """
    Parse MHC alleles into their canonical representation.

    Examples:
        A2 -> HLA-A2
        A0201 -> HLA-A*02:01
        H2-K-k -> H2-Kk
        RT-1*9.5:f -> RT1-9.5f
        DRB1_0101 -> HLA-DRB1*01:01

    Parameters
    ----------
    raw_string : str
        String corresponding to allele, locus, or other MHC-related name

    include_species_prefix : bool
        Include species in the normalized. If False, then you would
        get "A*02:01" for "A0201", instead of "HLA-A*02:01"

    infer_class2_pairing : bool
        If given only the alpha or beta chain of a Class II allele,
        try to infer the most likely pairing from population frequencies.

    default_species_prefix : str
        By default, parse alleles like "A*02:01" as human but it's possible
        to change this to some other species.
    """
    parsed_object = parse(
        raw_string,
        infer_class2_pairing=infer_class2_pairing,
        default_species_prefix=default_species_prefix)
    result = parsed_object.normalized_string(
        include_species=include_species_prefix)
    print(">>>", raw_string, parsed_object, result)
    return result


def compact_string(
        raw_string,
        infer_class2_pairing=False,
        default_species_prefix="HLA"):

    """
    Turn HLA-A*02:01 into A0201 or H-2-D-b into H-2Db or
    HLA-DPA1*01:05-DPB1*100:01 into DPA10105-DPB110001

    Parameters
    ----------
    raw_string : str
        String corresponding to allele, locus, or other MHC-related name

    infer_class2_pairing : bool
        If given only the alpha or beta chain of a Class II allele,
        try to infer the most likely pairing from population frequencies.

    default_species_prefix : str
        By default, parse alleles like "A*02:01" as human but it's possible
        to change this to some other species.
    """
    parsed_object = parse(
        raw_string,
        infer_class2_pairing=infer_class2_pairing,
        default_species_prefix=default_species_prefix)
    return parsed_object.compact_string(include_species_prefix=True)


def dataframe_from_list(names, default_species_prefix="HLA"):
    parsed_objects = [
        parse(name, default_species_prefix=default_species_prefix)
        for name in names
    ]
    records = [
        obj.to_dict()
        for obj in parsed_objects
    ]
    return pd.DataFrame.from_records(records)
