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

from collections import defaultdict
from .helpers import normalize_string


class NormalizingDictionary(object):
    """
    Like a regular dictionary but all keys get normalized by a user
    provided function.

    Caution: the number of items in keys() and values() for this dictionary
    may differ because two distinct keys may be transformed to the same
    underlying normalized key and thus will share a value.
    """
    def __init__(
            self,
            *pairs,
            normalize_fn=normalize_string,
            default_value_fn=None):
        self.store = {}
        self.original_to_normalized_key_dict = {}
        self.normalized_to_original_keys_dict = defaultdict(set)
        self.normalize_fn = normalize_fn
        self.default_value_fn = default_value_fn

        # populate dictionary with initial values via calls to __setitem__
        for (k, v) in pairs:
            self[k] = v

    def __getitem__(self, k):
        k_normalized = self.normalize_fn(k)
        if k_normalized not in self.store:
            if self.default_value_fn is not None:
                self.store[k_normalized] = self.default_value_fn()
            else:
                raise KeyError(k)
        return self.store[k_normalized]

    def __setitem__(self, k, v):
        k_normalized = self.normalize_fn(k)
        self.original_to_normalized_key_dict[k] = k_normalized
        self.normalized_to_original_keys_dict[k_normalized].add(k)
        self.store[k_normalized] = v

    def __contains__(self, k):
        k_normalized = self.normalize_fn(k)
        return k_normalized in self.store

    def original_keys(self, k):
        """
        Returns set of original keys which match the normalized representation
        of the given key.
        """
        k_normalized = self.normalize_fn(k)
        return self.normalized_to_original_keys_dict.get(k_normalized, set())

    def original_key(self, k):
        """
        Attempts to get a single original key matching which has
        the same the normalized representation as the given input,
        but may raise an exception if there are more than one original
        key that matches.
        """
        ks = list(self.original_keys(k))
        if len(ks) == 0:
            raise KeyError(k)
        elif len(ks) > 1:
            raise ValueError("Key '%s' matches multiple entries: %s" % (ks,))
        else:
            return ks[0]

    def get(self, k, v=None):
        return self.store.get(self.normalize_fn(k), v)

    def keys(self):
        return self.original_to_normalized_keys_dict.keys()

    def normalized_keys(self):
        return self.store.keys()

    def keys_aligned_with_values(self):
        """
        Returns one of the original keys associated with each item
        of values
        """
        return [
            self.normalized_to_original_keys_dict[k]
            for k in self.normalized_keys()
        ]

    def values(self):
        return self.store.values()

    def items(self):
        return zip(
            self.keys_aligned_with_values(),
            self.values())

    def map_values(self, fn):
        pairs = [
            (k, fn(v))
            for (k, v) in self.items()
        ]
        return NormalizingDictionary(
            *pairs,
            normalize_fn=self.normalize_fn,
            default_value_fn=self.default_value_fn)

    def map_keys(self, fn):
        pairs = [
            (fn(k), v)
            for (k, v) in self.items()
        ]
        return NormalizingDictionary(
            *pairs,
            normalize_fn=self.normalize_fn,
            default_value_fn=self.default_value_fn)

    def invert(self):
        """
        Returns a NormalizingDictionary where every value
        is associated with a set of keys which mapped to it. When
        values are collections (list, set, or tuple) then the elements
        of the collection are turned into individual keys.
        """
        result = NormalizingDictionary(
            normalize_fn=self.normalize_fn,
            default_value_fn=set)

        for (k, v) in self.items():
            if isinstance(v, (list, tuple, set)):
                values = v
            else:
                values = [v]
            for vi in values:
                result[vi].add(k)
        return result

    @classmethod
    def from_dict(cls, d, normalize_fn=normalize_string, default_value_fn=None):
        return cls(
            *d.items(),
            normalize_fn=normalize_fn,
            default_value_fn=default_value_fn)

    def __str__(self):
        s = (
            "<NormalizingDictionary with %d unique items, "
            "normalize_fn=%s, "
            "default_value_fn=%s>") % (
                len(self.store),
                self.normalize_fn,
                self.default_value_fn)
        all_items = list(self.items())

        for i, (k, v) in enumerate(all_items):
            s += "\n\t%s: %s" % (k, v)
            if i > 10:
                s += "\n..."
                break
        return s

    def __repr__(self):
        return str(self)
