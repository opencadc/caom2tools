=========
Caom2Repo
=========

.. image:: https://img.shields.io/pypi/v/caom2repo.svg
    :target: https://pypi.python.org/pypi/caom2repo

Client caom2-repo
=================
caom2Repo provides a client (caom2-repo) to perform CRUD (Create, Read, Update, Delete) on an observation in a collection in a repository.

Visitor Plugin
==============
The client also provides a visitor function which accepts a plugin. The visitor function iterates the observations of a collection and updates them according to the algorithm of the plugin function. The following is an example plugin to add a 'PREVIEW' Plane to an observation. More plugin examples can be found in caom2repo/tests/. ::

    from __future__ import (absolute_import, division, print_function,
                            unicode_literals)

    from caom2.observation import Observation
    from caom2.plane import Plane


    class ObservationUpdater(object):
        """ObservationUpdater that adds a plane to the observation."""

        def update(self, observation, **kwargs):
            """
            Processes an observation and updates it
            """
            assert isinstance(observation, Observation), (
                "observation %s is not an Observation".format(observation))
            observation.planes.add(Plane('PREVIEW'))
