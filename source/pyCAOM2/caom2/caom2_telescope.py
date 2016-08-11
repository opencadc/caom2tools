#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#***********************************************************************
#******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
#*************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2010.                            (c) 2010.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
#***********************************************************************
#

"""Definition of the caom2.telescope class"""

import util.caom2_util  as util
from caom2_object import Caom2Object


class Telescope(Caom2Object):
    """ Telescope """

    def __init__(self, name,
                 geo_location_x=None,
                 geo_location_y=None,
                 geo_location_z=None,
                 keywords=None
                 ):
        """
        Initializes a Telescope instance

        Arguments:
        name : name of the telescope
        """

        assert name is not None, "No telescope name provided"
        assert isinstance(name, str), "name is not a str: {0}".format(name)
        self.name = name
        self.geo_location_x = geo_location_x
        self.geo_location_y = geo_location_y
        self.geo_location_z = geo_location_z
        if keywords is None:
            keywords = set()
        self.keywords = keywords

    # Properties

    @property
    def name(self):
        """a name for this facility.

        eg. CFHT
        type: str

        """
        return self._name

    @name.setter
    def name(self, value):
        util.typeCheck(value, str, 'name', override=False)
        self._name = value

    @property
    def keywords(self):
        """A set that contains keywords associated with this telescope

        eg.  keywords.add('big')
        type: set

        """
        return self._keywords

    @keywords.setter
    def keywords(self, value):
        util.typeCheck(value, set, 'keywords', override=False)
        self._keywords = value

    @property
    def geo_location_x(self):
        """The x geocentric location of the telescope.

        This should be valid at time of MJD-OBS.

        These coordinates should be in the ITRS reference,
        basically whatever a GPS device is saying.

        The directions of the x/y/z follow the
        Earth Centred Rotation frame.

        units: m
        type: float

        """
        return self._geo_location_x

    @geo_location_x.setter
    def geo_location_x(self, value):
        util.typeCheck(value, float, 'geo_location_x')
        self._geo_location_x = value

    @property
    def geo_location_y(self):
        """the y geocentric (ECR) location of the telescope.

        This should be valid at time of MJD-OBS.

        These coordinates should be in the ITRS reference,
        basically whatever a GPS device is saying.

        The directions of the x/y/z follow the
        Earth Centred Rotation frame.

        units: m
        type: float

        """
        return self._geo_location_y

    @geo_location_y.setter
    def geo_location_y(self, value):
        util.typeCheck(value, float, 'geo_location_y')
        self._geo_location_y = value

    @property
    def geo_location_z(self):
        """the z geocentric (ECR) location of the telescope.
        This should be valid at time of MJD-OBS.

        These coordinates should be in the ITRS reference,
        basically whatever a GPS device is saying.

        The directions of the x/y/z follow the
        Earth Centred Rotation frame.

        units: m
        type: float

        """
        return self._geo_location_z

    @geo_location_z.setter
    def geo_location_z(self, value):
        util.typeCheck(value, float, 'geo_location_z')
        self._geo_location_z = value
