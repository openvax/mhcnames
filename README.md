> **⚠️ DEPRECATED:** mhcnames is no longer maintained. Please use [mhcgnomes](https://github.com/openvax/mhcgnomes) instead, which provides a superset of mhcnames functionality with support for more MHC nomenclature systems.

---

# **DEPRECATED: you should use [mhcgnomes](https://pypi.org/project/mhcgnomes/) instead**

[![Tests](https://github.com/openvax/mhcnames/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/mhcnames/actions/workflows/tests.yml)
<a href="https://pypi.python.org/pypi/mhcnames/">
    <img src="https://img.shields.io/pypi/v/mhcnames.svg?maxAge=1000" alt="PyPI" />
</a>

# mhcnames
All the fun and adventure of MHC nomenclature, now in Python

Python 3.11 or later is required. New projects should use `mhcgnomes` instead.

## Example

```python
In [1]: mhcnames.parse_allele_name("HLA-A0201")
Out[1]: AlleleName(species='HLA', gene='A', allele_family='02', allele_code='01')

In [2]: mhcnames.compact_allele_name("HLA-A*02:01")
Out[2]: 'A0201'
```
