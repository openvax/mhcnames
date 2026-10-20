# Releasing mhcnames

This package is deprecated in favor of [mhcgnomes](https://github.com/openvax/mhcgnomes).
Compatibility releases still use the following checks:

1. Bump `mhcnames.__version__` in a pull request.
2. Run `./lint.sh` and `./test.sh`, then merge only after the Python 3.11–3.14 CI matrix passes.
3. From clean, up-to-date `master`, build and check both distributions:

   ```sh
   python -m build
   python -m twine check --strict dist/*
   ```

4. Upload the distributions to PyPI, tag the merged commit with the version, and verify the published metadata.
