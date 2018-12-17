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


def expand_strings(*names, chars_to_remove="-_'"):
    """
    Returns set of names which includes original input, lowercase, uppercase,
    and case variants of strings with dash, underline, and apostrophe removed.
    """
    result_set = set(names)

    for char in chars_to_remove:
        for name in list(result_set):
            result_set.add(name.replace(char, ""))

    for name in list(result_set):
        result_set.add(name.upper())

    for name in list(result_set):
        result_set.add(name.lower())
    return result_set


def apply_string_expansion_to_set_members(set_of_strings):
    """
    For every string in the given set, include all of its uppercase and
    dash-less variants in the result set.
    """
    result = set([])
    for x in set_of_strings:
        for y in expand_strings(x):
            result.add(y)
    return result


def apply_string_expansion_to_dict_keys(d):
    """
    Create a larger dictionary by copying value associated with each key
    to uppercase and dash-less variants of the key.
    """
    result = {}
    for key, value in d.items():
        for expanded_key in expand_strings(key):
            result[expanded_key] = value
    return result
