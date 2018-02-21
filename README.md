<a href="https://travis-ci.org/openvax/mhcnames">
    <img src="https://travis-ci.org/openvax/mhcnames.svg?branch=master" alt="Build Status" />
</a>
<a href="https://coveralls.io/github/openvax/mhcnames?branch=master">
    <img src="https://coveralls.io/repos/openvax/mhcnames/badge.svg?branch=master&service=github" alt="Coverage Status" />
</a>
<a href="https://pypi.python.org/pypi/mhcnames/">
    <img src="https://img.shields.io/pypi/v/mhcnames.svg?maxAge=1000" alt="PyPI" />
</a>

# mhcnames
All the fun and adventure of MHC nomenclature, now in Python

## Example

```python
In [1]: mhcnames.parse_allele_name("HLA-A0201")
Out[1]: AlleleName(species='HLA', gene='A', allele_family='02', allele_code='01')

In [2]: mhcnames.compact_allele_name("HLA-A*02:01")
Out[2]: 'A0201'
```
