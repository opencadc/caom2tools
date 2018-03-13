caom2utils
=========

.. image:: https://img.shields.io/pypi/v/caom2utils.svg   
    :target: https://pypi.python.org/pypi/caom2utils

Utilities to facililate working with the CAOM2 model.

Observation Validation
---------------------

Validates a CAOM2 element (Observation, Plane, Artifact, Part or Chunk) with respect to the attributes of the element and possibly all its sub-elements. Example of validations: attribute values, spherical geometry of planes, WCS of chuncks, etc.

.. code:: python

    from __future__ import (absolute_import, division, print_function,
                    unicode_literals)
    import sys
    from caom2 import SimpleObservation
    import caom2utils

    obs = SimpleObservation('collection', 'observationID')

    # change and update obs

    try:
        caom2utils.validate(obs)
    except Exception:
        print('My exception is not valid')
