#
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

""" Deines __init__ """

#
from caom_object import CaomObject
from caom_object import AbstractCaomEntity

# Util classes
from util import Util
from util import TypedList
from util import TypedSet
from util import TypedOrderedDict
from util import ClassProperty
from util import Validator

# WCS data types
from wcs import Axis
from wcs import Coord2D
from wcs import CoordAxis1D
from wcs import CoordAxis2D
from wcs import CoordBounds1D
from wcs import CoordBounds2D
from wcs import CoordCircle2D
from wcs import CoordError
from wcs import CoordFunction1D
from wcs import CoordFunction2D
from wcs import CoordPolygon2D
from wcs import CoordRange1D
from wcs import CoordRange2D
from wcs import Dimension2D
from wcs import RefCoord
from wcs import Slice
from wcs import ValueCoord2D
from wcs import EnergyTransition

# Discovery data types
from shape import Box
from shape import Circle
from shape import Interval
from shape import Point
from shape import Polygon
from shape import Vertex

# Chunk level classes
from chunk import Chunk
from chunk import ObservableAxis
from chunk import SpatialWCS
from chunk import SpectralWCS
from chunk import TemporalWCS
from chunk import PolarizationWCS

# Part level classes
from part import Part

# Artifact level classes
from artifact import Artifact

# Plane level classes
from plane import Plane
from plane import PlaneURI
from plane import DataQuality
from plane import Metrics
from plane import Provenance
from plane import Position
from plane import Energy
from plane import EnergyTransition
from plane import Polarization
from plane import Time

# Observation level classes
from observation import Observation
from observation import ObservationURI
from observation import SimpleObservation
from observation import CompositeObservation
from observation import Algorithm
from observation import Environment
from observation import Proposal
from observation import Requirements
from observation import Target
from observation import TargetPosition
from observation import Telescope

# enums
from artifact import ProductType
from artifact import ReleaseType
from plane import CalibrationLevel
from plane import DataProductType
from plane import EnergyBand
from plane import PolarizationState
from plane import Quality
from observation import ObservationIntentType
from observation import Status
from observation import TargetType

# observation reader and writer
from obs_reader_writer import ObservationReader
from obs_reader_writer import ObservationWriter
from obs_reader_writer import CAOM20_NAMESPACE
from obs_reader_writer import CAOM21_NAMESPACE
from obs_reader_writer import CAOM22_NAMESPACE
