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


def parse_substring(allele, pred, max_len=None):
    """
    Extract substring of letters for which predicate is True
    """
    result = ""
    pos = 0
    if max_len is None:
        max_len = len(allele)
    else:
        max_len = min(max_len, len(allele))
    while pos < max_len and pred(allele[pos]):
        result += allele[pos]
        pos += 1
    return result, allele[pos:]


SEPARATORS = {":", "*", "-"}

def parse_separator(allele, max_len=None):
    return parse_substring(allele, lambda c: c in SEPARATORS, max_len=max_len)

def parse_alphanum(allele, max_len=None):
    return parse_substring(allele, lambda c: c.isalnum(), max_len=max_len)

def parse_letters(allele, max_len=None):
    return parse_substring(allele, lambda c: c.isalpha(), max_len=max_len)

def parse_numbers(allele, max_len=None):
    return parse_substring(allele, lambda c: c.isdigit(), max_len=max_len)

def parse_not_numbers(allele, max_len=None):
    return parse_substring(allele, lambda c: not c.isdigit(), max_len=max_len)

def parse_until(allele, sep):
    return parse_substring(allele, lambda c: c != sep)
