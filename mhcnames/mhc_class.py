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

class1_subtypes = {
    "Ia",
    "Ib",
    "Ic"
    "Id"
}

class2_subtypes = {
    "IIa",
    "IIb",
    "IIc"
    "IId"
}

valid_class_restrictions = {
    "I",
    "II",
}.union(class1_subtypes).union(class2_subtypes)

classical_subtypes = {"Ia", "IIa"}


def mhc_class_is_more_specific(original_class=None, new_class=None):
    return (
        (original_class is None and new_class is not None) or
        (original_class == "I" and new_class in class1_subtypes) or
        (original_class == "II" and new_class in class1_subtypes)
    )


def is_valid_restriction(original_class=None, new_class=None):
    if original_class is None:
        return True

    if new_class is None:
        # once we've restricted, can't go backwards
        return False

    if original_class == new_class:
        return True

    return mhc_class_is_more_specific(original_class, new_class)


def restrict_alleles(alleles, mhc_class):
    if mhc_class == "I":
        valid_subtypes = class1_subtypes
    elif mhc_class == "II":
        valid_subtypes = class2_subtypes
    else:
        valid_subtypes = {mhc_class}

    return [
        allele
        for allele in alleles
        if allele.get_mhc_class() in valid_subtypes
    ]
