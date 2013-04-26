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

""" defines the time class"""


from types.caom2_interval import Interval
import util.caom2_util as util

class Time(object):
    """ Time """

    def __init__(self,
                 value=None,
                 bounds=None,
                 dimension=None,
                 resolution=None,
                 sample_size=None,
                 exposure=None):
        """
        Initialize a Time instance.

        Arguments:
        None
        """
        self.value = value
        self.bounds = bounds
        self.dimension = dimension
        self.resolution = resolution
        self.sample_size = sample_size
        self.exposure = exposure

    # Properties

    @property
    def value(self):
        """ Actual time value, seconds since epoch.

        quesion is, what epoch?

        units: s
        """
        return self._value

    @value.setter
    def value(self, value):
        util.typeCheck(value, float, 'value')
        self._value = value

    @property
    def bounds(self):
        """an interval object that gives start and end of an time interval

        Not actually implemented yet.  Likely you want a TemporalWCS
        instead

        type: Interval(lower_mjd, upper_mjd)
        unit: mjd
        
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        util.typeCheck(value, Interval, 'bounds')
        self._bounds = value

    @property
    def dimension(self):
        """Number of pixel in the time direction, normally 1.
        
        eg 1
        type: long

        """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        util.typeCheck(value, long, 'dimension')
        self._dimension = value

    @property
    def resolution(self):
        """Time resolution of the samples, in seconds.

        normally this is the same as the exposure time,
        but in a stack you might have a larger resolution value than
        exposure time.

        eg. 1000
        unit: s
        type: float

        """
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        util.typeCheck(value, float, 'resolution')
        self._resolution = value

    @property
    def sample_size(self):
        """nominally the exposure time, in seconds.

        """
        return self._sample_size

    @sample_size.setter
    def sample_size(self, value):
        util.typeCheck(value, float, 'sample_size')
        self._sample_size = value

    @property
    def exposure(self):
        """Duration of the exposure, in seconds"""
        return self._exposure

    @exposure.setter
    def exposure(self, value):
        util.typeCheck(value, float, 'exposure')
        self._exposure = value
