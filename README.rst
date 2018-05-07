Bio2BEL CTD
===========
Enrich BEL graphs with the impact of chemical perturbagens on biological entities and systems.

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
``bio2bel_ctd`` can be installed easily from `PyPI <https://pypi.python.org/pypi/bio2bel_ctd>`_ with the
following code in your favorite terminal:

.. code-block:: sh

    $ python3 -m pip install bio2bel_ctd

or from the latest code on `GitHub <https://github.com/bio2bel/ctd>`_ with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/bio2bel/ctd.git@master

Setup
-----
The CTD can be downloaded and populated from either the Python REPL or the automatically installed command line
utility.

Python REPL
~~~~~~~~~~~
.. code-block:: python

    >>> import bio2bel_ctd
    >>> ctd_manager = bio2bel_ctd.Manager()
    >>> ctd_manager.populate()

Command Line Utility
~~~~~~~~~~~~~~~~~~~~
.. code-block:: bash

    bio2bel_ctd populate

Acknowledgements
----------------
This package heavily relies on Christian Ebeling's `PyCTD <https://github.com/cebel/pyctd>`_.

.. |build| image:: https://travis-ci.org/bio2bel/ctd.svg?branch=master
    :target: https://travis-ci.org/bio2bel/ctd
    :alt: Build Status

.. |coverage| image:: https://codecov.io/gh/bio2bel/ctd/coverage.svg?branch=master
    :target: https://codecov.io/gh/bio2bel/ctd?branch=master
    :alt: Coverage Status

.. |documentation| image:: https://readthedocs.org/projects/ctd/badge/?version=latest
    :target: http://ctd.readthedocs.io
    :alt: Documentation Status

.. |climate| image:: https://codeclimate.com/github/bio2bel/ctd/badges/gpa.svg
    :target: https://codeclimate.com/github/bio2bel/ctd
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/bio2bel_ctd.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/bio2bel_ctd.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/bio2bel_ctd.svg
    :alt: MIT License
