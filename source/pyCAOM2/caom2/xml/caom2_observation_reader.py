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

""" Defines ObservationReader class """

from lxml import etree
import pkg_resources
from .. caom2_algorithm import Algorithm
from .. caom2_artifact import Artifact
from .. caom2_enums import DataProductType
from .. caom2_enums import CalibrationLevel
from .. caom2_enums import ProductType
from .. caom2_enums import ObservationIntentType
from .. caom2_enums import TargetType
from .. caom2_chunk import Chunk
from .. caom2_composite_observation import CompositeObservation
from .. caom2_energy_transition import EnergyTransition
from .. caom2_environment import Environment
from .. caom2_exceptions import ObservationParsingException
from .. caom2_instrument import Instrument
from .. caom2_metrics import Metrics
from .. caom2_observation_uri import ObservationURI
from .. caom2_part import Part
from .. caom2_plane import Plane
from .. caom2_plane_uri import PlaneURI
from .. caom2_proposal import Proposal
from .. caom2_provenance import Provenance
from .. caom2_simple_observation import SimpleObservation
from .. caom2_target import Target
from .. caom2_target_position import TargetPosition
from .. caom2_telescope import Telescope
from .. util.caom2_util import str2ivoa
from .. wcs.caom2_axis import Axis
from .. wcs.caom2_coord2d import Coord2D
from .. wcs.caom2_value_coord2d import ValueCoord2D
from .. wcs.caom2_coord_axis1d import CoordAxis1D
from .. wcs.caom2_coord_axis2d import CoordAxis2D
from .. wcs.caom2_coord_bounds1d import CoordBounds1D
from .. wcs.caom2_coord_circle2d import CoordCircle2D
from .. wcs.caom2_coord_error import CoordError
from .. wcs.caom2_coord_function1d import CoordFunction1D
from .. wcs.caom2_coord_function2d import CoordFunction2D
from .. wcs.caom2_coord_polygon2d import CoordPolygon2D
from .. wcs.caom2_coord_range1d import CoordRange1D
from .. wcs.caom2_coord_range2d import CoordRange2D
from .. wcs.caom2_dimension2d import Dimension2D
from .. wcs.caom2_observable_axis import ObservableAxis
from .. wcs.caom2_polarization_wcs import PolarizationWCS
from .. wcs.caom2_ref_coord import RefCoord
from .. wcs.caom2_slice import Slice
from .. wcs.caom2_spatial_wcs import SpatialWCS
from .. wcs.caom2_spectral_wcs import SpectralWCS
from .. wcs.caom2_temporal_wcs import TemporalWCS
from .. types.caom2_point import Point


class ObservationReader(object):
    """ObservationReader """

    CAOM2_PKG = 'caom2'
    SCHEMA_FILE = 'CAOM-2.0.xsd'

    def __init__(self, valididate):
        """Constructor. XML Schema validation may be disabled, in which case
        the client is likely to fail in horrible ways if it received invalid
        documents. However, performance may be improved.

        Arguments:
        validate : True if enable schema validation, False otherwise
        """
        self._valididate = valididate
        schema_path = pkg_resources.resource_filename(
            ObservationReader.CAOM2_PKG, ObservationReader.SCHEMA_FILE)
        xmlschema_doc = etree.parse(schema_path)
        self._xmlschema = etree.XMLSchema(xmlschema_doc)

    def _getChildElement(self, elTag, parent, ns, required):
        for element in list(parent):
            if (element.tag == "{" + ns + "}" + elTag):
                if (not element.keys() and not element.text):
                    # element is empty, return None
                    return None
                else:
                    # element has content, return it
                    return element

        if (required):
            error = elTag + " element not found in " + parent.tag
            raise ObservationParsingException(error)
        else:
            return None

    def _getChildText(self, elTag, parent, ns, required):
        childElement = self._getChildElement(elTag, parent, ns, required)
        if (childElement is None):
            return None
        else:
            return childElement.text

    def _getChildTextAsInt(self, elTag, parent, ns, required):
        childElement = self._getChildElement(elTag, parent, ns, required)
        if (childElement is None):
            return None
        else:
            return int(childElement.text)

    def _getChildTextAsLong(self, elTag, parent, ns, required):
        childElement = self._getChildElement(elTag, parent, ns, required)
        if (childElement is None):
            return None
        else:
            return long(childElement.text)

    def _getChildTextAsFloat(self, elTag, parent, ns, required):
        childElement = self._getChildElement(elTag, parent, ns, required)
        if (childElement is None):
            return None
        else:
            return float(childElement.text)

    def _getAlgorithm(self, elTag, parent, ns, required):
        """Build an Algorithm object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the algorithm element
        ns : namespace of the document
        required : indicates whether the element is required
        return : an Algorithm object or
                 None if the document does not contain an algorithm element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return Algorithm(self._getChildText("name", el, ns, True))

    def _getMetaRelease(self, elTag, parent, ns, required):
        """Build a MetaRelease object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the MetaRelease element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a MetaRelease object or
                None if the document does not contain a MetaRelease element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            # TODO: need to catch exceptions,
            # what kind of exceptions are thrown?
            return str2ivoa(el.text)

    def _getProposal(self, elTag, parent, ns, required):
        """Build a Proposal object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Proposal element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Proposal object or
                 None if the document does not contain a Proposal element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            proposal = Proposal(self._getChildText("id", el, ns, True))
            proposal.pi_name = self._getChildText("pi", el, ns, False)
            proposal.project = self._getChildText("project", el, ns, False)
            proposal.title = self._getChildText("title", el, ns, False)
            keywords = self._getChildText("keywords", el, ns, False)
            if (keywords is not  None):
                proposal.keywords.list = keywords.split()
            return proposal

    def _getTarget(self, elTag, parent, ns, required):
        """Build a Target object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Target element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Target object or
                 None if the document does not contain a Target element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            target = Target(self._getChildText("name", el, ns, True))
            targetType = self._getChildText("type", el, ns, False)
            if (targetType):
                target.target_type = TargetType.getByValue(targetType)
            target.standard = ("true" ==
                self._getChildText("standard", el, ns, False))
            target.redshift = (
                self._getChildTextAsFloat("redshift", el, ns, False))
            target.moving = ("true" ==
                self._getChildText("moving", el, ns, False))
            keywords = self._getChildText("keywords", el, ns, False)
            if (keywords is not None):
                target.keywords.list = keywords.split()
            return target

    def _getTargetPosition(self, elTag, parent, ns, required):
        """Build a TargetPosition object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Target element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a TargetPosition object or
                 None if the document does not contain a TargetPosition element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            target_position = TargetPosition(
                self._getPoint("coordinates", el, ns, True),
                self._getChildText("coordsys", el, ns, True))
            target_position.equinox = (
                self._getChildTextAsFloat("equinox", el, ns, False))
            return target_position

    def _getTelescope(self, elTag, parent, ns, required):
        """Build a Telescope object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Telescope element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Telescope object or
                 None if the document does not contain a Telescope element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            telescope = Telescope(self._getChildText("name", el, ns, True))
            telescope.geo_location_x = (
                self._getChildTextAsFloat("geoLocationX", el, ns, False))
            telescope.geo_location_y = (
                self._getChildTextAsFloat("geoLocationY", el, ns, False))
            telescope.geo_location_z = (
                self._getChildTextAsFloat("geoLocationZ", el, ns, False))
            keywords = self._getChildText("keywords", el, ns, False)
            if (keywords is not None):
                telescope.keywords.list = keywords.split()
            return telescope

    def _getInstrument(self, elTag, parent, ns, required):
        """Build an Instrument object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Instrument element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Instrument object or
                 None if the document does not contain an Instrument element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            instrument = Instrument(self._getChildText("name", el, ns, True))
            keywords = self._getChildText("keywords", el, ns, False)
            if (keywords is not None):
                instrument.keywords.list = keywords.split()
            return instrument

    def _getEnvironment(self, elTag, parent, ns, required):
        """Build an Environment object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Environment element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Environment object or
                 None if the document does not contain an Environment element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            environment = Environment()
            environment.seeing = (
                self._getChildTextAsFloat("seeing", el, ns, False))
            environment.humidity = (
                self._getChildTextAsFloat("humidity", el, ns, False))
            environment.elevation = (
                self._getChildTextAsFloat("elevation", el, ns, False))
            environment.tau = (
                self._getChildTextAsFloat("tau", el, ns, False))
            environment.wavelength_tau = (
                self._getChildTextAsFloat("wavelengthTau", el, ns, False))
            environment.ambient_temp = (
                self._getChildTextAsFloat("ambientTemp", el, ns, False))
            environment.photometric = ("true" ==
                self._getChildText("photometric", el, ns, False))
            return environment

    def _addMembers(self, members, parent, ns):
        """Create ObservationURI objects from an XML representation of
        ObservationURI elements found in members element, and add them to the
        set of ObservationURI's

        Arguments:
        members : Set of member's from the parent Observation object
        parent : element containing the Environment element
        ns : namespace of the document
        return : an Environment object or
                 None if the document does not contain an Environment element
        raise : ObservationParsingException
        """
        el = self._getChildElement("members", parent, ns, False)
        if (el is not None):
            for memberEl in el.iterchildren("{" + ns + "}observationURI"):
                members.add(ObservationURI(memberEl.text))

            if not members:
                error = "No observationURI element found in members"
                raise ObservationParsingException(error)

    def _addInputs(self, inputs, parent, ns):
        """Create PlaneURI objects from an XML representation of the planeURI
        elements and add them to the set of PlaneURIs.

        Arguments:
        inputs : set of PlaneURI from the Provenance
        parent : element containing the PlaneURI elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._getChildElement("inputs", parent, ns, False)
        if (el is not None):
            for uriEl in el.iterchildren("{" + ns + "}planeURI"):
                inputs.add(PlaneURI(uriEl.text))

            if not inputs:
                error = "No planeURI element found in members"
                raise ObservationParsingException(error)

    def _getProvenance(self, elTag, parent, ns, required):
        """Build a Provenance object from an XML representation of a
        Provenance element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Provenance element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Provenance object or
                 None if the document does not contain a Provenance element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            prov = Provenance(self._getChildText("name", el, ns, True))
            prov.version = self._getChildText("version", el, ns, False)
            prov.project = self._getChildText("project", el, ns, False)
            prov.producer = self._getChildText("producer", el, ns, False)
            prov.run_id = self._getChildText("runID", el, ns, False)
            reference = self._getChildText("reference", el, ns, False)
            if reference:
                prov.reference = reference
            prov.last_executed = str2ivoa(
                self._getChildText("lastExecuted", el, ns, False))
            keywords = self._getChildText("keywords", el, ns, False)
            if (keywords is not None):
                prov.keywords.list = keywords.split()
            self._addInputs(prov.inputs, el, ns)
            return prov

    def _getMetrics(self, elTag, parent, ns, required):
        """Build a Metrics object from an XML representation of a
        Metrics element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Metrics element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Metrics object or
                 None if the document does not contain a Metrics element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            metrics = Metrics()
            metrics.source_number_density = \
                self._getChildTextAsFloat("sourceNumberDensity", el, ns, False)
            metrics.background = \
                self._getChildTextAsFloat("background", el, ns, False)
            metrics.background_std_dev = \
                self._getChildTextAsFloat("backgroundStddev", el, ns, False)
            metrics.flux_density_limit = \
                self._getChildTextAsFloat("fluxDensityLimit", el, ns, False)
            metrics.mag_limit = \
                self._getChildTextAsFloat("magLimit", el, ns, False)
            return metrics

    def _getPoint(self, elTag, parent, ns, required):
        """Build an Point object from an XML representation
        of an Point element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Point element
        ns : namespace of the document
        required : indicate whether the element is required
        return : an Point object or
                 None if the document does not contain an Point element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return Point(self._getChildTextAsFloat("cval1", el, ns, True),
                         self._getChildTextAsFloat("cval2", el, ns, True))

    def _getAxis(self, elTag, parent, ns, required):
        """Build an Axis object from an XML representation of an Axis element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Axis element
        ns : namespace of the document
        required : indicate whether the element is required
        return : an Axis object or
                 None if the document does not contain an Axis element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return Axis(self._getChildText("ctype", el, ns, True),
                        self._getChildText("cunit", el, ns, False))

    def _getSlice(self, elTag, parent, ns, required):
        """Build a Slice object from an XML representation of a Slice element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Slice element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a Slice object or
                 None if the document does not contain a Slice element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return Slice(self._getAxis("axis", el, ns, True),
                         self._getChildTextAsLong("bin", el, ns, True))

    def _getObservableAxis(self, elTag, parent, ns, required):
        """Build an ObservableAxis object from an XML representation of an
        observable element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Observable element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : an ObservableAxis object or
                 None if the document does not contain an Observable element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            observable = ObservableAxis(
                self._getSlice("dependent", el, ns, True))
            observable.independent = \
                self._getSlice("independent", el, ns, False)
            return observable

    def _getCoordError(self, elTag, parent, ns, required):
        """Build a CoordError object from an XML representation of an error
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the error element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordError object or
                 None if the document does not contain an error element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return CoordError(
                self._getChildTextAsFloat("syser", el, ns, True),
                self._getChildTextAsFloat("rnder", el, ns, True))

    def _getRefCoord(self, elTag, parent, ns, required):
        """Build a RefCoord object from an XML representation of a coord
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the coord element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a RefCoord object or
                 None if the document does not contain a coord element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return RefCoord(
                self._getChildTextAsFloat("pix", el, ns, True),
                self._getChildTextAsFloat("val", el, ns, True))

    def _getCoord2D(self, elTag, parent, ns, required):
        """Build a Coord2D object from an XML representation of a coord
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the coord element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a Coord2D object or
                 None if the document does not contain a coord element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return Coord2D(
                self._getRefCoord("coord1", el, ns, True),
                self._getRefCoord("coord2", el, ns, True))

    def _getValueCoord2D(self, elTag, parent, ns, required):
        """Build a ValueCoord2D object from an XML representation of a
        value coord element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the coord element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a ValueCoord2D object or
                 None if the document does not contain a coord element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return ValueCoord2D(
                self._getChildTextAsFloat("coord1", el, ns, True),
                self._getChildTextAsFloat("coord2", el, ns, True))

    def _getCoordRange2D(self, elTag, parent, ns, required):
        """Build a CoordRange2D object from an XML representation of a range
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the range element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordRange2D object or
                 None if the document does not contain a range element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return CoordRange2D(
                self._getCoord2D("start", el, ns, True),
                self._getCoord2D("end", el, ns, True))

    def _getCoordCircle2D(self, elTag, parent, ns, required):
        """Build a CoordCircle2D object from an XML representation of a circle
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the circle element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordCircle2D object or
                 None if the document does not contain a circle element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return CoordCircle2D(
                self._getValueCoord2D("center", el, ns, True),
                self._getChildTextAsFloat("radius", el, ns, True))

    def _getCoordPolygon2D(self, elTag, parent, ns, required):
        """Build a CoordPolygon2D object from an XML representation
        of a polygon element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the polygon element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordPolygon2D object or
                 None if the document does not contain a polygon element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            verticeEl = self._getChildElement("vertices", el, ns, True)
            childrenVertices = list(
                verticeEl.iterchildren(tag=("{" + ns + "}vertex")))
            if (len(childrenVertices) < 3):
                error = ("CoordPolygon2D must have a minimum of 3 vertices, "
                    "found " + len(childrenVertices))
                raise ObservationParsingException(error)
            else:
                polygon = CoordPolygon2D()
                for childVertexEl in childrenVertices:
                    polygon.vertices.append(ValueCoord2D(
                        self._getChildTextAsFloat(
                                    "coord1", childVertexEl, ns, True),
                        self._getChildTextAsFloat(
                                    "coord2", childVertexEl, ns, True)))
                return polygon

    def _getCoordBounds2D(self, elTag, parent, ns, required):
        """Build a CoordBounds2D object from an XML representation of a bounds
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the bounds element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordBounds2D object or
                 None if the document does not contain a bounds element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            circle = self._getCoordCircle2D("circle", el, ns, False)
            if (circle is not None):
                return circle
            else:
                polygon = self._getCoordPolygon2D("polygon", el, ns, False)
                if (polygon is not None):
                    return polygon
                else:
                    error = "Unsupported element not found in " + elTag + \
                        ": " + el.getText()
                    raise ObservationParsingException(error)

    def _getDimension2D(self, elTag, parent, ns, required):
        """Build a Dimension2D object from an XML representation of a dimension
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the dimension element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a Dimention2D object or
                 None if the document does not contain a dimension element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return Dimension2D(
                self._getChildTextAsLong("naxis1", el, ns, True),
                self._getChildTextAsLong("naxis2", el, ns, True))

    def _getCoordFunction2D(self, elTag, parent, ns, required):
        """Build a CoordFunction2D object from an XML representation of a
        function element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the function element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordFunction2D object or
                 None if the document does not contain a function element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return CoordFunction2D(
                self._getDimension2D("dimension", el, ns, True),
                self._getCoord2D("refCoord", el, ns, True),
                self._getChildTextAsFloat("cd11", el, ns, True),
                self._getChildTextAsFloat("cd12", el, ns, True),
                self._getChildTextAsFloat("cd21", el, ns, True),
                self._getChildTextAsFloat("cd22", el, ns, True))

    def _getCoordAxis2D(self, elTag, parent, ns, required):
        """Build a CoordAxis2D object from an XML representation of an axis
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the axis element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordAxis2D object or
                 None if the document does not contain an axis element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            axis = CoordAxis2D(self._getAxis("axis1", el, ns, True),
                               self._getAxis("axis2", el, ns, True))
            axis.error1 = self._getCoordError("error1", el, ns, False)
            axis.error2 = self._getCoordError("error2", el, ns, False)
            axis.range = self._getCoordRange2D("range", el, ns, False)
            axis.bounds = self._getCoordBounds2D("bounds", el, ns, False)
            axis.function = self._getCoordFunction2D("function", el, ns, False)
            return axis

    def _getSpatialWCS(self, elTag, parent, ns, required):
        """Build a SpatialWCS object from an XML representation of a position
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a SpatialWCS object or
                 None if the document does not contain a position element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            position = SpatialWCS(self._getCoordAxis2D("axis", el, ns, False))
            position.coordsys = self._getChildText("coordsys", el, ns, False)
            position.equinox = \
                self._getChildTextAsFloat("equinox", el, ns, False)
            position.resolution = \
                self._getChildTextAsFloat("resolution", el, ns, False)
            return position

    def _addChildrenToCoordRange1DList(self, elTag, ranges, parent, ns,
                                       required):
        """Create CoordRange1D objects from an XML representation of the
        range elements and add them to the set of ranges.

        Arguments:
        elTag : element tag which identifies the element
        ranges : reference to set of ranges
        parent : element containing the ranges elements
        ns : namespace of the document
        required : boolean indicating whether the element is required
        """
        for rangeEl in parent.iterchildren("{" + ns + "}" + elTag):
            ranges.append(CoordRange1D(
                 self._getRefCoord("start", rangeEl, ns, True),
                 self._getRefCoord("end", rangeEl, ns, True)))

    def _getCoordBounds1D(self, elTag, parent, ns, required):
        """Build a CoordBounds1D object from an XML representation of a bounds
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the bounds element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordBounds1D object or
                 None if the document does not contain a bounds element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            coordBounds1D = CoordBounds1D()
            samplesEl = self._getChildElement("samples", el, ns, False)
            if (samplesEl is not None):
                self._addChildrenToCoordRange1DList(
                    "range", coordBounds1D.samples, samplesEl, ns, False)
            return coordBounds1D

    def _getCoordRange1D(self, elTag, parent, ns, required):
        """Build a CoordRange1D object from an XML representation of a range
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the range element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordRange1D object or
                 None if the document does not contain a range element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return CoordRange1D(
                self._getRefCoord("start", el, ns, True),
                self._getRefCoord("end", el, ns, True))

    def _getCoordFunction1D(self, elTag, parent, ns, required):
        """Build a CoordFunction1D object from an XML representation of a
        function element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the function element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordFunction1D object or
                 None if the document does not contain a function element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return CoordFunction1D(
                self._getChildTextAsLong("naxis", el, ns, True),
                self._getChildTextAsFloat("delta", el, ns, True),
                self._getRefCoord("refCoord", el, ns, True))

    def _getCoordAxis1D(self, elTag, parent, ns, required):
        """Build a CoordAxis1D object from an XML representation of an axis
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the axis element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CoordAxis1D object or
                 None if the document does not contain an axis element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            axis = CoordAxis1D(self._getAxis("axis", el, ns, True))
            axis.error = self._getCoordError("error", el, ns, False)
            axis.range = self._getCoordRange1D("range", el, ns, False)
            axis.bounds = self._getCoordBounds1D("bounds", el, ns, False)
            axis.function = self._getCoordFunction1D("function", el, ns, False)
            return axis

    def _getTransition(self, elTag, parent, ns, required):
        """Build an EnergyTransition object from an XML representation of
        a transition element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the transition element
        ns : namespace of the document
        required : boolean indicating whether the element is reuiqred
        return : an EnergyTransition object or
                 None if the document does not contain a transition element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return EnergyTransition(
                self._getChildText("species", el, ns, True),
                self._getChildText("transition", el, ns, True))

    def _getSpectralWCS(self, elTag, parent, ns, required):
        """Build a SpectralWCS object from an XML representation of an energy
        element.

        Arguments:
        elTag : element tag which indentifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a SpectralWCS object or
                 None if the document does not contain an energy element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            energy = SpectralWCS(
                self._getCoordAxis1D("axis", el, ns, True),
                self._getChildText("specsys", el, ns, True))
            energy.ssysobs = \
                self._getChildText("ssysobs", el, ns, False)
            energy.ssyssrc = \
                self._getChildText("ssyssrc", el, ns, False)
            energy.restfrq = \
                self._getChildTextAsFloat("restfrq", el, ns, False)
            energy.restwav = \
                self._getChildTextAsFloat("restwav", el, ns, False)
            energy.velosys = \
                self._getChildTextAsFloat("velosys", el, ns, False)
            energy.zsource = \
                self._getChildTextAsFloat("zsource", el, ns, False)
            energy.velang = \
                self._getChildTextAsFloat("velang", el, ns, False)
            energy.bandpassName = \
                self._getChildText("bandpassName", el, ns, False)
            energy.resolvingPower = \
                self._getChildTextAsFloat("resolvingPower", el, ns, False)
            energy.transition = \
                self._getTransition("transition", el, ns, False)
            return energy

    def _getTemporalWCS(self, elTag, parent, ns, required):
        """Build a TemporalWCS object from an XML representation of an time
        element.

        Arguments:
        elTag : element tag which indentifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a TemporalWCS object or
                 None if the document does not contain an time element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            time = TemporalWCS(
                self._getCoordAxis1D("axis", el, ns, True))
            time.timesys = \
                self._getChildText("timesys", el, ns, False)
            time.trefpos = \
                self._getChildText("trefpos", el, ns, False)
            time.mjdref = \
                self._getChildTextAsFloat("mjdref", el, ns, False)
            time.exposure = \
                self._getChildTextAsFloat("exposure", el, ns, False)
            time.resolution = \
                self._getChildTextAsFloat("resolution", el, ns, False)
            return time

    def _getPolarizationWCS(self, elTag, parent, ns, required):
        """Build a PolarizationWCS object from an XML representation of a
        polarization element.

        Arguments:
        elTag : element tag which indentifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a PolarizationWCS object or
                 None if the document does not contain an polarization element
        raise : ObservationParsingException
        """
        el = self._getChildElement(elTag, parent, ns, required)
        if (el is None):
            return None
        else:
            return PolarizationWCS(
                self._getCoordAxis1D("axis", el, ns, False))

    def _addChunks(self, chunks, parent, ns):
        """Build Chunk objects from an XML representation of Chunk elements
        and add them to the set of Chunks.

        Argument:
        chunks : set of Chunk objects from the Part
        parent : element containing the Chunk elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._getChildElement("chunks", parent, ns, False)
        if (el is None):
            return None
        else:
            for chunkEl in el.iterchildren("{" + ns + "}chunk"):
                tempChunk = Chunk()
                tempChunk._id = \
                    long(chunkEl.get("{" + ns + "}id"))
                tempChunk._last_modified = \
                    str2ivoa(chunkEl.get(
                        "{" + ns + "}lastModified"))
                productType = \
                    self._getChildText("productType", chunkEl, ns, False)
                if (productType):
                    tempChunk.product_type = \
                        ProductType.getByValue(productType)
                tempChunk.naxis = \
                    self._getChildTextAsInt("naxis", chunkEl, ns, False)
                tempChunk.observable_axis = \
                    self._getChildTextAsInt("observableAxis", chunkEl, ns,
                                            False)
                tempChunk.position_axis_1 = \
                    self._getChildTextAsInt("positionAxis1", chunkEl, ns,
                                            False)
                tempChunk.position_axis_2 = \
                    self._getChildTextAsInt("positionAxis2", chunkEl, ns,
                                            False)
                tempChunk.energy_axis = \
                    self._getChildTextAsInt("energyAxis", chunkEl, ns, False)
                tempChunk.time_axis = \
                    self._getChildTextAsInt("timeAxis", chunkEl, ns, False)
                tempChunk.polarization_axis = \
                    self._getChildTextAsInt("polarizationAxis", chunkEl, ns,
                                            False)
                tempChunk.observable = \
                    self._getObservableAxis("observable", chunkEl, ns, False)
                tempChunk.position = \
                    self._getSpatialWCS("position", chunkEl, ns, False)
                tempChunk.energy = \
                    self._getSpectralWCS("energy", chunkEl, ns, False)
                tempChunk.time = \
                    self._getTemporalWCS("time", chunkEl, ns, False)
                tempChunk.polarization = \
                    self._getPolarizationWCS("polarization", chunkEl, ns,
                                             False)
                chunks.append(tempChunk)

    def _addParts(self, parts, parent, ns):
        """Build Part objects from an XML representation of Part elements and
        add them to the set of Parts.

        Argument:
        parts : set of Part objects from the Artifact
        parent : element containing the Part elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._getChildElement("parts", parent, ns, False)
        if (el is None):
            return None
        else:
            for partEl in el.iterchildren("{" + ns + "}part"):
                tempPart = \
                    Part(self._getChildText("name", partEl, ns, True))
                tempPart._id = \
                    long(partEl.get("{" + ns + "}id"))
                tempPart._last_modified = \
                    str2ivoa(partEl.get(
                        "{" + ns + "}lastModified"))
                productType = \
                    self._getChildText("productType", partEl, ns, False)
                if (productType):
                    tempPart.product_type = \
                        ProductType.getByValue(productType)
                self._addChunks(tempPart.chunks, partEl, ns)
                parts[tempPart.name] = tempPart

    def _addArtifacts(self, artifacts, parent, ns):
        """Build artifacts from an XML representation of the artifact elements
        and add them to the set of Artifacts.

        Arguments:
        artifacts : set of Artifacts from the Plane
        parent : element containing the Artifact elements
        ns : namespace fo the document
        raise : ObservationParsingException
        """
        el = self._getChildElement("artifacts", parent, ns, False)
        if (el is None):
            return None
        else:
            for artifactEl in el.iterchildren("{" + ns + "}artifact"):
                tempArtifact = \
                    Artifact(self._getChildText("uri", artifactEl, ns, True))
                tempArtifact._id = \
                    long(artifactEl.get("{" + ns + "}id"))
                tempArtifact._last_modified = \
                    str2ivoa(artifactEl.get(
                        "{" + ns + "}lastModified"))
                tempArtifact.content_type = \
                    self._getChildText("contentType", artifactEl, ns, False)
                tempArtifact.content_length = \
                    self._getChildTextAsLong("contentLength", artifactEl, ns,
                                             False)
                productType = \
                    self._getChildText("productType", artifactEl, ns, False)
                if (productType):
                    tempArtifact.product_type = \
                        ProductType.getByValue(productType)
                tempArtifact.alternative = "true" == (
                    self._getChildText("alternative", artifactEl, ns, False))
                self._addParts(tempArtifact.parts, artifactEl, ns)
                artifacts[tempArtifact.uri] = tempArtifact

    def _addPlanes(self, planes, parent, ns):
        """Create Planes object from XML representation of Plane elements
        and add them to the set of Planes.

        Arguments:
        planes : Set of planes from the parent Observation object
        parent : element containing the Plane elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._getChildElement("planes", parent, ns, False)
        if (el is None):
            return None
        else:
            for planeEl in el.iterchildren("{" + ns + "}plane"):
                tempPlane = Plane(
                    self._getChildText("productID", planeEl, ns, True))
                tempPlane._id = \
                    long(planeEl.get("{" + ns + "}id"))
                tempPlane._last_modified = \
                    str2ivoa(planeEl.get(
                        "{" + ns + "}lastModified"))
                tempPlane.meta_release = str2ivoa(
                    self._getChildText("metaRelease", planeEl, ns, False))
                tempPlane.data_release = str2ivoa(
                    self._getChildText("dataRelease", planeEl, ns, False))
                dataProductType = \
                    self._getChildText("dataProductType", planeEl, ns, False)
                if (dataProductType):
                    tempPlane.data_product_type = \
                        DataProductType.getByValue(dataProductType)
                calibrationLevel = \
                    self._getChildText("calibrationLevel", planeEl, ns, False)
                if (calibrationLevel):
                    tempPlane.calibration_level = \
                        CalibrationLevel.getByValue(int(calibrationLevel))
                tempPlane.provenance = \
                    self._getProvenance("provenance", planeEl, ns, False)
                tempPlane.metrics = \
                    self._getMetrics("metrics", planeEl, ns, False)
                self._addArtifacts(tempPlane.artifacts, planeEl, ns)
                planes[tempPlane.product_id] = tempPlane

            if not planes:
                error = "No plane element found in planes"
                raise ObservationParsingException(error)

    def read(self, source):
        """Build an Observation object from an XML document located in source.
        Source an be a file name/path, a file object, a file-like object or a
        URL using the HTTP or FTP protocol.

        Arguments:
        source : source of XML document containing an Observation element
        return : an Observation object
        raise : ObservationParsingException
        """
        doc = etree.parse(source)
        if self._valididate:
            self._xmlschema.validate(doc)
        self._xmlschema.assert_(etree.parse(source))
        root = doc.getroot()
        ns = root.nsmap["caom2"]
        collection = self._getChildElement("collection", root, ns, True).text
        observationID = \
            self._getChildElement("observationID", root, ns, True).text
        # Instantiate Algorithm
        algorithm = self._getAlgorithm("algorithm", root, ns, True)
        # Instantiate Observation
        if (root.get("{http://www.w3.org/2001/XMLSchema-instance}type")
            == "caom2:SimpleObservation"):
            observation = SimpleObservation(collection, observationID)
            observation.algorithm = algorithm
        else:
            observation = \
                CompositeObservation(collection, observationID, algorithm)
        # Instantiate children of Observation
        observation._id = long(root.get("{" + ns + "}id"))
        observation._last_modified = \
            str2ivoa(root.get("{" + ns + "}lastModified"))
        observation.sequence_number = \
            self._getChildTextAsInt("sequenceNumber", root, ns, False)
        observation.obs_type = \
            self._getChildText("type", root, ns, False)
        intent = self._getChildText("intent", root, ns, False)
        if (intent):
            observation.intent = ObservationIntentType.getByValue(intent)
        observation.meta_release = \
            self._getMetaRelease("metaRelease", root, ns, False)
        observation.proposal = \
            self._getProposal("proposal", root, ns, False)
        observation.target = \
            self._getTarget("target", root, ns, False)
        observation.target_position = \
            self._getTargetPosition("targetPosition", root, ns, False)
        observation.telescope = \
            self._getTelescope("telescope", root, ns, False)
        observation.instrument = \
            self._getInstrument("instrument", root, ns, False)
        observation.environment = \
            self._getEnvironment("environment", root, ns, False)
        self._addPlanes(observation.planes, root, ns)
        if isinstance(observation, CompositeObservation):
            self._addMembers(observation.members, root, ns)

        return observation

if __name__ == '__main__':
    reader = ObservationReader(True)
    reader.read('../../../../test/data/CompleteCompositePolygon.xml')
