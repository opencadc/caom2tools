Common Archive Observation Model - data engineering tools

.. image:: https://img.shields.io/pypi/pyversions/caom2.svg
    :target: https://pypi.python.org/pypi/caom2

.. image:: https://github.com/opencadc/caom2tools/workflows/CI/badge.svg?branch=master&event=schedule
    :target: https://github.com/opencadc/caom2tools/actions?query=event%3Aschedule+

.. image:: https://codecov.io/gh/opencadc/caom2tools/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/opencadc/caom2tools

.. image:: https://img.shields.io/github/contributors/opencadc/caom2tools.svg
    :target: https://github.com/opencadc/caom2tools/graphs/contributors

Set of Python tools for working with the cADC CAOM2 data model: https://www.opencadc.org/caom2/


Developers Guide
================


Requires pip.

Installing Packages
-------------------
Note: might need to escape chars in your shell

::

    cd caom2 && pip install -e .[test]
    cd caom2utils && pip install -e .[test]
    cd caom2repo && pip install -e .[test]

Testing packages
----------------

Testing caom2
~~~~~~~~~~~~~~~~~

::

    cd ./caom2
    pytest caom2

Testing caom2utils
~~~~~~~~~~~~~~~~

::

    cd ./caom2utils
    pytest caom2utils

Testing caom2rep
~~~~~~~~~~~~~~~~

::

    cd ./caom2repo
    pytest caom2repo



Checkstyle
~~~~~~~~~~
flake8 style checking is enforced on pull requests. Following commands should
not report errors

::

     flake8 caom2/caom2 caom2utils/caom2utils caom2repo/caom2repo


Testing with tox
~~~~~~~~~~~~~~~~

If tox, the generic virtual environment tool, is available it can be used to test with different versions of
python is isolation. For example, to test on all supported versions of Python in cadcdata (assuming that
they are available in the system):

::

    cd ./caom2repo && tox

To test a specific version:

::

    cd ./caom2utils && tox -e py3.9


To list all the available environments:

::

    cd ./caom2 && tox -a

