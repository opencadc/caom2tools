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

"""definition of the caom2.Metrics object"""

from caom2_object import Caom2Object
import util.caom2_util as util


class Metrics(Caom2Object):
    """ Metrics """

    def __init__(self):
        """
        Initializes a Metrics instance

        Arguments:
        None
        """
        self._source_number_density = None
        self._background = None
        self._background_std_dev = None
        self._flux_density_limit = None
        self._mag_limit = None

    # Properties
    @property
    def source_number_density(self):
        """The number of sources brighter than mag_limit (flux_density_limit)

        unit: ct/deg2
        type: float
        """
        return self._source_number_density

    @source_number_density.setter
    def source_number_density(self, value):
        util.typeCheck(value, float, "source_number_density")
        util.valueCheck(value, 0, 1E10, "source_number_density")
        self._source_number_density = value

    @property
    def background(self):
        """The flux in the sky (background).

        units: Jy/pix
        type: float
        """
        return self._background

    @background.setter
    def background(self, value):
        util.typeCheck(value, float, "background")
        util.valueCheck(value, 0, 1E10, "background")
        self._background = value

    @property
    def background_std_dev(self):
        """the standard deviation (per pixel) in background flux.

        Likely this only makes sense to define if background is also defined.
        units: Jy/pix
        type: float
        """
        return self._background_std_dev

    @background_std_dev.setter
    def background_std_dev(self, value):
        util.typeCheck(value, float, "background_std_dev")
        util.valueCheck(value, 0, 1E10, "background")
        self._background_std_dev = value

    @property
    def flux_density_limit(self):
        """flux density where S:N=5 for point source.

        this is intended to provide a measure of the limit of detection.

        units: Jy
        type: float
        """
        return self._flux_density_limit

    @flux_density_limit.setter
    def flux_density_limit(self, value):
        util.typeCheck(value, float, "flux_denisty_limit")
        util.valueCheck(value, 0, 1E10, "flux_density_limit")
        self._flux_density_limit = value

    @property
    def mag_limit(self):
        """AB magnitude limit where S:N=5 for point source.

        Likely specify just mag_limit or flux_density_limit, not both?

        units: AB mag
        type: float
        """
        return self._mag_limit

    @mag_limit.setter
    def mag_limit(self, value):
        util.typeCheck(value, float, 'mag_limit')
        util.valueCheck(value, 0, 40, 'mag_limit')
        self._mag_limit = value
