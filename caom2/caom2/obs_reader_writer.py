# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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
# ***********************************************************************
#


""" Defines ObservationReader class """
import json
import os
import uuid
from builtins import str, int
from enum import Enum
from io import IOBase

from lxml import etree

import caom2
from . import artifact
from . import caom_util
from . import chunk
from . import dali
from . import observation
from . import part
from . import plane
from . import shape
from . import wcs
import logging

from .plane import CalibrationStatus, Ucd, Polarization

DATA_PKG = 'data'

CAOM23_SCHEMA_FILE = 'CAOM-2.3.xsd'
CAOM24_SCHEMA_FILE = 'CAOM-2.4.xsd'
CAOM25_SCHEMA_FILE = 'CAOM-2.5.xsd'

CAOM23_NAMESPACE = 'http://www.opencadc.org/caom2/xml/v2.3'
CAOM24_NAMESPACE = 'http://www.opencadc.org/caom2/xml/v2.4'
CAOM25_NAMESPACE = 'http://www.opencadc.org/caom2/xml/v2.5'

CAOM_VERSION = {
                CAOM23_NAMESPACE: 23,
                CAOM24_NAMESPACE: 24,
                CAOM25_NAMESPACE: 25
                }

CAOM23 = "{%s}" % CAOM23_NAMESPACE
CAOM24 = "{%s}" % CAOM24_NAMESPACE
CAOM25 = "{%s}" % CAOM25_NAMESPACE

XSI_NAMESPACE = "http://www.w3.org/2001/XMLSchema-instance"
XSI = "{%s}" % XSI_NAMESPACE

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

__all__ = ['ObservationReader', 'ObservationWriter',
           'ObservationParsingException']

logger = logging.getLogger(__name__)


def _to_samples(vertices):
    samples = []
    last_closed_point = None
    points = []
    for vertex in vertices:
        if vertex.type == shape.SegmentType.MOVE:
            points.append(shape.Point(vertex.cval1, vertex.cval2))
        elif vertex.type == shape.SegmentType.CLOSE:
            last_closed_point = shape.Point(vertex.cval1, vertex.cval2)
            points.append(last_closed_point)
            samples.append(shape.Polygon(points))
            points = []  # continue with a new polygon
        else:
            if not points:
                # no move so start from the last closed point
                points.append(last_closed_point)
            points.append(shape.Point(vertex.cval1, vertex.cval2))
    return caom2.MultiShape(samples)


def _to_vertices(samples):
    vertices = []
    for sample_shape in samples.shapes:
        if isinstance(sample_shape, shape.Polygon):
            vertices.append(shape.Vertex(sample_shape.points[0].cval1,
                                         sample_shape.points[0].cval2,
                                         shape.SegmentType.MOVE))
            for point in sample_shape.points[1:-1]:
                vertices.append(shape.Vertex(point.cval1, point.cval2,
                                             shape.SegmentType.LINE))
            vertices.append(shape.Vertex(sample_shape.points[-1].cval1,
                                         sample_shape.points[-1].cval2,
                                         shape.SegmentType.CLOSE))
        else:
            raise ValueError("Only polygons can be converted to Vertices/MultiPolygon")

    return vertices


class ObservationReader(object):
    """ObservationReader """

    def __init__(self, validate=False, input_format='xml'):
        """Constructor. XML Schema validation may be disabled, in which case
        the client is likely to fail in horrible ways if it received invalid
        documents. However, performance may be improved.

        Arguments:
        validate : If True enable schema validation, False otherwise
        """
        if input_format not in ['xml', 'json']:
            raise ValueError("Unsupported format: " + input_format)
        self._input_format = input_format
        self._validate = validate
        if self._validate:
            if input_format == 'json':
                raise ValueError("Schema validation not supported for JSON")
            caom23_schema_path = os.path.join(THIS_DIR + '/' + DATA_PKG,
                                              CAOM23_SCHEMA_FILE)

            parser = etree.XMLParser(remove_blank_text=True)
            xsd = etree.parse(caom23_schema_path, parser)

            caom24_schema = etree.Element(
                '{http://www.w3.org/2001/XMLSchema}import',
                namespace=CAOM24_NAMESPACE,
                schemaLocation=CAOM24_SCHEMA_FILE)
            xsd.getroot().insert(1, caom24_schema)

            caom25_schema = etree.Element(
                '{http://www.w3.org/2001/XMLSchema}import',
                namespace=CAOM25_NAMESPACE,
                schemaLocation=CAOM25_SCHEMA_FILE)
            xsd.getroot().insert(2, caom25_schema)

            self._xmlschema = etree.XMLSchema(xsd)
            self.version = None

    def _get_attribute(self, name, ns, element):
        if self._input_format == 'xml':
            return element.get("{" + ns + "}" + name)
        else:
            return element['@' + name]

    def _set_entity_attributes(self, element, ns, caom2_entity):
        element_id = self._get_attribute("id", ns, element)
        element_last_modified = self._get_attribute("lastModified", ns, element)
        element_max_last_modified = self._get_attribute("maxLastModified", ns, element)
        element_meta_checksum = self._get_attribute("metaChecksum", ns, element)
        element_acc_meta_checksum = self._get_attribute("accMetaChecksum", ns, element)
        element_meta_producer = self._get_attribute("metaProducer", ns, element)

        caom2_entity._id = uuid.UUID(element_id)

        if element_last_modified:
            caom2_entity._last_modified = caom_util.str2ivoa(
                element_last_modified)
        if element_max_last_modified:
            caom2_entity._max_last_modified = caom_util.str2ivoa(
                element_max_last_modified)
        if element_meta_checksum:
            caom2_entity._meta_checksum = element_meta_checksum
        if element_acc_meta_checksum:
            caom2_entity._acc_meta_checksum = element_acc_meta_checksum
        if element_meta_producer:
            caom2_entity._meta_producer = element_meta_producer

    def _get_child_element(self, element_tag, parent, ns, required):
        if not parent:
            if required:
                error = "Parent element is None, cannot find " + element_tag
                raise ObservationParsingException(error)
            return None
        if self._input_format == 'xml':
            for element in list(parent):
                if element.tag == "{" + ns + "}" + element_tag:
                    if list(element) == 0 and not element.keys() and\
                       (not element.text or not element.text.strip()):
                        # element is empty, return None
                        return None
                    else:
                        # element has content, return it
                        return element
            if required:
                error = element_tag + " element not found in " + parent.tag
                raise ObservationParsingException(error)
        elif element_tag in parent:
            if element_tag in parent:
                return parent[element_tag]
            if required:
                error = element_tag + " element not found in " + parent.uri
                raise ObservationParsingException(error)

            return None

    def _get_child_text(self, element_tag, parent, ns, required):
        child_element = self._get_child_element(element_tag, parent, ns,
                                                required)
        if child_element is None:
            return None
        if self._input_format == 'xml':
            return str(child_element.text)
        else:
            return str(child_element)

    def _get_child_text_as_int(self, element_tag, parent, ns, required):
        child_element = self._get_child_element(element_tag, parent, ns,
                                                required)
        if child_element is None:
            return None

        if self._input_format == 'xml':
            return int(child_element.text)
        else:
            return int(child_element)

    def _get_child_text_as_long(self, element_tag, parent, ns, required):
        return self._get_child_text_as_int(element_tag, parent, ns, required)

    def _get_child_text_as_float(self, element_tag, parent, ns, required):
        child_element = self._get_child_element(element_tag, parent, ns,
                                                required)
        if child_element is None:
            return None
        if self._input_format == 'xml':
            return float(child_element.text)
        else:
            return float(child_element)

    def _get_child_text_as_boolean(self, element_tag, parent, ns, required):
        child_element = self._get_child_element(element_tag, parent, ns,
                                                required)
        if child_element is None:
            return None
        if self._input_format == 'xml':
            return child_element.text.lower() == "true"
        else:
            return child_element

    def _add_keywords(self, keywords_list, element, ns, required):
        """
        Parses keywords sub-elements and adds them to a list

        :param keywords_list: list to add the keywords to
        :param element: xml element
        :param ns: name space
        :param required: keywords sub-element required or not
        """
        keywords_element = self._get_child_element("keywords", element, ns,
                                                   required)
        if self._input_format == 'xml':
            if keywords_element is not None:
                for keyword in self._get_children(keywords_element, ns, "keyword"):
                    keywords_list.add(keyword.text)
        else:
            for keyword in keywords_element:
                keywords_list.add(keyword)

    def _get_algorithm(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return observation.Algorithm(
                self._get_child_text("name", el, ns, True))

    def _get_meta_release(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            # TODO: need to catch exceptions,
            # what kind of exceptions are thrown?

            if self._input_format == 'xml':
                return caom_util.str2ivoa(el.text)
            else:
                return caom_util.str2ivoa(el)

    def _get_groups(self, element_tag, parent, ns, required):
        """Build set of groups from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the MetaRelease element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a URISet of group URIs or None
                None if the document does not contain meta groups element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            result = caom_util.URISet()
            if self._input_format == 'xml':
                for member_element in self._get_children(el, ns, "groupURI"):
                    result.add(member_element.text)
            else:
                for member in el:
                    result.add(member)
            return result

    def _get_proposal(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            proposal = observation.Proposal(
                self._get_child_text("id", el, ns, True))
            proposal.pi = self._get_child_text("pi", el, ns, False)
            proposal.project = self._get_child_text("project", el, ns, False)
            proposal.title = self._get_child_text("title", el, ns, False)
            self._add_keywords(proposal.keywords, el, ns, False)
            proposal.reference = self._get_child_text("reference", el, ns, False)
            return proposal

    def _get_target(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            target = observation.Target(
                self._get_child_text("name", el, ns, True))
            target_type = self._get_child_text("type", el, ns, False)
            if target_type:
                target.type = observation.TargetType(target_type)
            target_standard = self._get_child_text("standard", el, ns, False)
            if target_standard is not None:
                target.standard = ("true" == target_standard)
            target.redshift = (
                self._get_child_text_as_float("redshift", el, ns, False))
            target_moving = self._get_child_text("moving", el, ns, False)
            if target_moving is not None:
                target.moving = ("true" == target_moving)
            self._add_keywords(target.keywords, el, ns, False)
            target_id = self._get_child_text("targetID", el, ns, False)
            if target_id is not None:
                target.target_id = target_id
            return target

    def _get_target_position(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            coords_element = self._get_child_element("coordinates", el, ns,
                                                     required)
            target_position = observation.TargetPosition(
                self._get_point(coords_element, ns, True),
                self._get_child_text("coordsys", el, ns, True))
            target_position.equinox = (
                self._get_child_text_as_float("equinox", el, ns, False))
            return target_position

    def _get_requirements(self, element_tag, parent, ns, required):
        """Build an Requirements object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Requirements element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Requirements object or
                 None if the document does not contain an Requirements element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            flag = self._get_child_text("flag", el, ns, True)
            requirements = observation.Requirements(observation.Status(flag))
            return requirements

    def _get_telescope(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            telescope = observation.Telescope(
                self._get_child_text("name", el, ns, True))
            telescope.geo_location_x = (
                self._get_child_text_as_float("geoLocationX", el, ns, False))
            telescope.geo_location_y = (
                self._get_child_text_as_float("geoLocationY", el, ns, False))
            telescope.geo_location_z = (
                self._get_child_text_as_float("geoLocationZ", el, ns, False))
            self._add_keywords(telescope.keywords, el, ns, False)
            telescope.tracking_mode = self._get_child_text("trackingMode", el, ns, False)
            return telescope

    def _get_instrument(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            instrument = observation.Instrument(
                self._get_child_text("name", el, ns, True))
            self._add_keywords(instrument.keywords, el, ns, False)
            return instrument

    def _get_environment(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            environment = observation.Environment()
            environment.seeing = (
                self._get_child_text_as_float("seeing", el, ns, False))
            environment.humidity = (
                self._get_child_text_as_float("humidity", el, ns, False))
            environment.elevation = (
                self._get_child_text_as_float("elevation", el, ns, False))
            environment.tau = (
                self._get_child_text_as_float("tau", el, ns, False))
            environment.wavelength_tau = (
                self._get_child_text_as_float("wavelengthTau", el, ns, False))
            environment.ambient_temp = (
                self._get_child_text_as_float("ambientTemp", el, ns, False))
            photometric = self._get_child_text("photometric", el, ns, False)
            if photometric is not None:
                environment.photometric = ("true" == photometric.lower())
            return environment

    def _add_members(self, members, parent, ns):
        """Create observation URI objects from an XML representation of
        observation URI elements found in members element, and add them to the
        set of observation URI's

        Arguments:
        members : Set of member's from the parent Observation object
        parent : element containing the Environment element
        ns : namespace of the document
        return : an Environment object or
                 None if the document does not contain an Environment element
        raise : ObservationParsingException
        """
        el = self._get_child_element("members", parent, ns, False)
        if el is not None:
            if self.version < 25:
                obs_tag = "observationURI"
            else:
                obs_tag = "member"
            for member_element in self._get_children(el, ns, obs_tag):
                if self._input_format == 'xml':
                    members.add(member_element.text)
                else:
                    members.add(member_element)

    def _add_inputs(self, inputs, parent, ns):
        """Create URI objects from an XML representation of the plane
        elements and add them to the set of plane URIs.

        Arguments:
        inputs : set of plane URIs from the Provenance
        parent : element containing the plane uri elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._get_child_element("inputs", parent, ns, False)
        if el is not None:
            if self.version < 25:
                plane_tag = "planeURI"
            else:
                plane_tag = "input"

            for uri_element in self._get_children(el, ns, plane_tag):
                if self._input_format == 'xml':
                    val = str(uri_element.text)
                else:
                    val = str(uri_element)
                inputs.add(val)

            if not inputs:
                error = "No planeURI element found in members"
                raise ObservationParsingException(error)

    def _get_provenance(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            prov = plane.Provenance(self._get_child_text("name", el, ns, True))
            prov.version = self._get_child_text("version", el, ns, False)
            prov.project = self._get_child_text("project", el, ns, False)
            prov.producer = self._get_child_text("producer", el, ns, False)
            prov.run_id = self._get_child_text("runID", el, ns, False)
            reference = self._get_child_text("reference", el, ns, False)
            if reference:
                prov.reference = reference
            prov.last_executed = caom_util.str2ivoa(
                self._get_child_text("lastExecuted", el, ns, False))
            self._add_keywords(prov.keywords, el, ns, False)
            self._add_inputs(prov.inputs, el, ns)
            return prov

    def _get_metrics(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            metrics = plane.Metrics()
            metrics.source_number_density = \
                self._get_child_text_as_float("sourceNumberDensity", el, ns,
                                              False)
            metrics.background = \
                self._get_child_text_as_float("background", el, ns, False)
            metrics.background_std_dev = \
                self._get_child_text_as_float("backgroundStddev", el, ns,
                                              False)
            metrics.flux_density_limit = \
                self._get_child_text_as_float("fluxDensityLimit", el, ns,
                                              False)
            metrics.mag_limit = \
                self._get_child_text_as_float("magLimit", el, ns, False)
            metrics.sample_snr = \
                self._get_child_text_as_float("sampleSNR", el, ns, False)
            return metrics

    def _get_quality(self, element_tag, parent, ns, required):
        """Build a Quality object from an XML representation

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the Quality element
        ns : namespace of the document
        required : indicates whether the element is required
        return : a Quality object or
                 None if the document does not contain an Quality element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            flag = self._get_child_text("flag", el, ns, True)
            data_quality = plane.DataQuality(plane.Quality(flag))
            return data_quality

    def _get_observable(self, parent, ns):
        """Build an Observable object from an XML representation

        Arguments:
        parent : element containing the Observable element
        ns : namespace of the document
        return : a Observable object or None if the document does not contain one
        raise : ObservationParsingException
        """
        el = self._get_child_element("observable", parent, ns, False)
        if el is None:
            return None
        else:
            ucd = self._get_child_text("ucd", el, ns, True)
            observable = plane.Observable(Ucd(ucd))
            calib = self._get_child_text("calibration", el, ns, False)
            if calib:
                observable.calibration = CalibrationStatus(calib)
            return observable

    def _get_point(self, point, ns, required):
        """Build an Point object from an XML representation
        of an Point element.

        Arguments:
        point : the point element
        ns : namespace of the document
        required : indicate whether the element is required
        return : an Point object
        raise : ObservationParsingException
        """
        return shape.Point(
            self._get_child_text_as_float("cval1", point, ns, True),
            self._get_child_text_as_float("cval2", point, ns, True))

    def _get_axis(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.Axis(self._get_child_text("ctype", el, ns, True),
                            self._get_child_text("cunit", el, ns, False))

    def _get_slice(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.Slice(self._get_axis("axis", el, ns, True),
                             self._get_child_text_as_long("bin", el, ns, True))

    def _get_observable_axis(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            observable = chunk.ObservableAxis(
                self._get_slice("dependent", el, ns, True))
            observable.independent = \
                self._get_slice("independent", el, ns, False)
            return observable

    def _get_coord_error(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.CoordError(
                self._get_child_text_as_float("syser", el, ns, True),
                self._get_child_text_as_float("rnder", el, ns, True))

    def _get_ref_coord(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.RefCoord(
                self._get_child_text_as_float("pix", el, ns, True),
                self._get_child_text_as_float("val", el, ns, True))

    def _get_coord2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.Coord2D(
                self._get_ref_coord("coord1", el, ns, True),
                self._get_ref_coord("coord2", el, ns, True))

    def _get_value_coord2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.ValueCoord2D(
                self._get_child_text_as_float("coord1", el, ns, True),
                self._get_child_text_as_float("coord2", el, ns, True))

    def _get_coord_range2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.CoordRange2D(
                self._get_coord2d("start", el, ns, True),
                self._get_coord2d("end", el, ns, True))

    def _get_coord_circle2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.CoordCircle2D(
                self._get_value_coord2d("center", el, ns, True),
                self._get_child_text_as_float("radius", el, ns, True))

    def _get_coord_polygon2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            vertice_element = self._get_child_element("vertices", el, ns, True)
            children_vertices = list(
                self._get_children(vertice_element, ns, "vertex"))
            if len(children_vertices) < 3:
                error = ("CoordPolygon2D must have a minimum of 3 vertices, "
                         "found " + len(children_vertices))
                raise ObservationParsingException(error)
            else:
                polygon = wcs.CoordPolygon2D()
                for child_vertex_el in children_vertices:
                    polygon.vertices.append(wcs.ValueCoord2D(
                        self._get_child_text_as_float(
                            "coord1", child_vertex_el, ns, True),
                        self._get_child_text_as_float(
                            "coord2", child_vertex_el, ns, True)))
                return polygon

    def _get_coord_bounds2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            circle = self._get_coord_circle2d("circle", el, ns, False)
            if circle is not None:
                return circle
            else:
                polygon = self._get_coord_polygon2d("polygon", el, ns, False)
                if polygon is not None:
                    return polygon
                else:
                    error = "Unsupported element not found in " + \
                            element_tag + ": " + el.getText()
                    raise ObservationParsingException(error)

    def _get_dimension2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.Dimension2D(
                self._get_child_text_as_long("naxis1", el, ns, True),
                self._get_child_text_as_long("naxis2", el, ns, True))

    def _get_coord_function2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.CoordFunction2D(
                self._get_dimension2d("dimension", el, ns, True),
                self._get_coord2d("refCoord", el, ns, True),
                self._get_child_text_as_float("cd11", el, ns, True),
                self._get_child_text_as_float("cd12", el, ns, True),
                self._get_child_text_as_float("cd21", el, ns, True),
                self._get_child_text_as_float("cd22", el, ns, True))

    def _get_coord_axis2d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            axis = wcs.CoordAxis2D(self._get_axis("axis1", el, ns, True),
                                   self._get_axis("axis2", el, ns, True))
            axis.error1 = self._get_coord_error("error1", el, ns, False)
            axis.error2 = self._get_coord_error("error2", el, ns, False)
            axis.range = self._get_coord_range2d("range", el, ns, False)
            axis.bounds = self._get_coord_bounds2d("bounds", el, ns, False)
            axis.function = self._get_coord_function2d("function", el, ns,
                                                       False)
            return axis

    def _get_spatial_wcs(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            position = chunk.SpatialWCS(
                self._get_coord_axis2d("axis", el, ns, False))
            position.coordsys = self._get_child_text("coordsys", el, ns, False)
            position.equinox = \
                self._get_child_text_as_float("equinox", el, ns, False)
            position.resolution = \
                self._get_child_text_as_float("resolution", el, ns, False)
            return position

    def _add_children_to_coord_range1d_list(self, element_tag, ranges, parent,
                                            ns, required):
        """Create CoordRange1D objects from an XML representation of the
        range elements and add them to the set of ranges.

        Arguments:
        elTag : element tag which identifies the element
        ranges : reference to set of ranges
        parent : element containing the ranges elements
        ns : namespace of the document
        required : boolean indicating whether the element is required
        """
        for range_element in self._get_children(parent, ns, element_tag):
            ranges.append(wcs.CoordRange1D(
                self._get_ref_coord("start", range_element, ns, True),
                self._get_ref_coord("end", range_element, ns, True)))

    def _get_coord_bounds1d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            coord_bounds1d = wcs.CoordBounds1D()
            samples_element = self._get_child_element("samples", el, ns, False)
            if samples_element is not None:
                self._add_children_to_coord_range1d_list(
                    "range", coord_bounds1d.samples, samples_element, ns,
                    False)
            return coord_bounds1d

    def _get_coord_range1d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.CoordRange1D(
                self._get_ref_coord("start", el, ns, True),
                self._get_ref_coord("end", el, ns, True))

    def _get_coord_function1d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.CoordFunction1D(
                self._get_child_text_as_long("naxis", el, ns, True),
                self._get_child_text_as_float("delta", el, ns, True),
                self._get_ref_coord("refCoord", el, ns, True))

    def _get_coord_axis1d(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            axis = wcs.CoordAxis1D(self._get_axis("axis", el, ns, True))
            axis.error = self._get_coord_error("error", el, ns, False)
            axis.range = self._get_coord_range1d("range", el, ns, False)
            axis.bounds = self._get_coord_bounds1d("bounds", el, ns, False)
            axis.function = self._get_coord_function1d("function", el, ns,
                                                       False)
            return axis

    def _get_transition(self, element_tag, parent, ns, required):
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
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return wcs.EnergyTransition(
                self._get_child_text("species", el, ns, True),
                self._get_child_text("transition", el, ns, True))

    def _get_spectral_wcs(self, element_tag, parent, ns, required):
        """Build a SpectralWCS object from an XML representation of an energy
        element.

        Arguments:
        elTag : element tag which indentifies the element
        parent : element containing the spectral wcs element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a SpectralWCS object or
                 None if the document does not contain an energy element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            energy = chunk.SpectralWCS(
                self._get_coord_axis1d("axis", el, ns, True),
                self._get_child_text("specsys", el, ns, True))
            energy.ssysobs = \
                self._get_child_text("ssysobs", el, ns, False)
            energy.ssyssrc = \
                self._get_child_text("ssyssrc", el, ns, False)
            energy.restfrq = \
                self._get_child_text_as_float("restfrq", el, ns, False)
            energy.restwav = \
                self._get_child_text_as_float("restwav", el, ns, False)
            energy.velosys = \
                self._get_child_text_as_float("velosys", el, ns, False)
            energy.zsource = \
                self._get_child_text_as_float("zsource", el, ns, False)
            energy.velang = \
                self._get_child_text_as_float("velang", el, ns, False)
            energy.bandpass_name = \
                self._get_child_text("bandpassName", el, ns, False)
            energy.resolving_power = \
                self._get_child_text_as_float("resolvingPower", el, ns, False)
            energy.transition = \
                self._get_transition("transition", el, ns, False)
            return energy

    def _get_temporal_wcs(self, element_tag, parent, ns, required):
        """Build a TemporalWCS object from an XML representation of an time
        element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the temporal wcs element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a TemporalWCS object or
                 None if the document does not contain an time element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            time = chunk.TemporalWCS(
                self._get_coord_axis1d("axis", el, ns, True))
            time.timesys = \
                self._get_child_text("timesys", el, ns, False)
            time.trefpos = \
                self._get_child_text("trefpos", el, ns, False)
            time.mjdref = \
                self._get_child_text_as_float("mjdref", el, ns, False)
            time.exposure = \
                self._get_child_text_as_float("exposure", el, ns, False)
            time.resolution = \
                self._get_child_text_as_float("resolution", el, ns, False)
            return time

    def _get_polarization_wcs(self, element_tag, parent, ns, required):
        """Build a PolarizationWCS object from an XML representation of a
        polarization element.

        Arguments:
        elTag : element tag which indentifies the element
        parent : element containing the polarization element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a PolarizationWCS object or
                 None if the document does not contain an polarization element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return chunk.PolarizationWCS(
                self._get_coord_axis1d("axis", el, ns, False))

    def _get_custom_wcs(self, element_tag, parent, ns, required):
        """Build a CustomWCS object from an XML representation of a
        polarization element.

        Arguments:
        elTag : element tag which indentifies the element
        parent : element containing the custom axis element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CustomWCS object or
                 None if the document does not contain a custom WCS element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        else:
            return chunk.CustomWCS(
                self._get_coord_axis1d("axis", el, ns, False))

    def _get_position(self, element_tag, parent, ns, required):
        """Build a Position object from an XML representation of a
        position element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a Position object or
                 None if the document does not contain a polarization element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        bounds, samples = self._get_shape("bounds", el, ns, False)
        pos = plane.Position(bounds=bounds, samples=samples)
        min_bounds = self._get_shape("minBounds", el, ns, False)
        if min_bounds:
            # ignore samples returned by get_shape
            pos.min_bounds = min_bounds[0]
        pos.dimension = self._get_dimension2d("dimension", el, ns, False)
        pos.max_recoverable_scale = self._get_interval("maxRecoverableScale", el, ns, False)
        pos.resolution = self._get_child_text_as_float("resolution", el, ns,
                                                       False)
        pos.resolution_bounds = self._get_interval("resolutionBounds", el,
                                                   ns, False)
        pos.sample_size = self._get_child_text_as_float("sampleSize", el, ns,
                                                        False)
        pos.calibration = self._get_child_text("calibration", el, ns, False)
        return pos

    def _get_energy(self, element_tag, parent, ns, required):
        """Build an Energy object from an XML representation of an
        energy element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : an Energy object or
                 None if the document does not contain an energy element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None

        bounds = self._get_interval("bounds", el, ns, True)
        if self.version < 25:
            samples = self._get_samples(self._get_child_element("bounds", el, ns, True), ns, True)
        else:
            samples = self._get_samples(el, ns, True)
        energy = plane.Energy(bounds, samples)
        energy.dimension = \
            self._get_child_text_as_int("dimension", el, ns, False)
        energy.resolving_power = self._get_child_text_as_float(
            "resolvingPower", el, ns, False)
        energy.resolving_power_bounds = self._get_interval(
            "resolvingPowerBounds", el, ns, False)
        energy.resolution = self._get_child_text_as_float("resolution", el, ns, False)
        energy.resolution_bounds = self._get_interval("resolutionBounds", el, ns, False)
        energy.sample_size = \
            self._get_child_text_as_float("sampleSize", el, ns, False)
        energy.bandpass_name = \
            self._get_child_text("bandpassName", el, ns, False)
        self._add_energy_bands(energy.energy_bands, el, ns)
        if self.version < 25:
            energy.rest = self._get_child_text_as_float("restwav", el, ns, False)
        else:
            energy.rest = \
                self._get_child_text_as_float("rest", el, ns, False)
        _transition_el = \
            self._get_child_element("transition", el, ns, required)
        if _transition_el is not None:
            species = self._get_child_text("species", _transition_el, ns, True)
            transition = \
                self._get_child_text("transition", _transition_el, ns, True)
            energy.transition = wcs.EnergyTransition(species, transition)
        energy.calibration = self._get_child_text("calibration", el, ns, False)

        return energy

    def _get_time(self, element_tag, parent, ns, required):
        """Build a Time object from an XML representation of a
        time element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a Time object or
                 None if the document does not contain a time element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        bounds = self._get_interval("bounds", el, ns, False)
        if self.version < 25:
            samples = self._get_samples(self._get_child_element("bounds", el, ns, True), ns, True)
        else:
            samples = self._get_samples(el, ns, True)
        time = plane.Time(bounds, samples)
        time.dimension = \
            self._get_child_text_as_int("dimension", el, ns, False)
        time.resolution = \
            self._get_child_text_as_float("resolution", el, ns, False)
        time.resolution_bounds = self._get_interval("resolutionBounds", el,
                                                    ns, False)
        time.sample_size = \
            self._get_child_text_as_float("sampleSize", el, ns, False)
        time.exposure = \
            self._get_child_text_as_float("exposure", el, ns, False)
        time.exposure_bounds = self._get_interval("exposureBounds", el, ns, False)
        time.calibration = self._get_child_text("calibration", el, ns, False)

        return time

    def _get_custom(self, element_tag, parent, ns, required):
        """Build a Custom object from an XML representation of a
        custom axis element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a CustomAxis object or
                 None if the document does not contain a custom axis element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        ctype = self._get_child_text("ctype", el, ns, False)
        bounds = self._get_interval("bounds", el, ns, False)
        if self.version < 25:
            samples = self._get_samples(self._get_child_element("bounds", el, ns, True), ns, True)
        else:
            samples = self._get_samples(el, ns, True)
        custom = plane.CustomAxis(ctype, bounds, samples=samples)
        custom.dimension = \
            self._get_child_text_as_int("dimension", el, ns, False)
        return custom

    def _get_visibility(self, parent, ns):
        """Build a Visibility object from an XML representation

        Arguments:
        parent : element containing the position element
        ns : namespace of the document
        return : a Visibility object or None if the document does not contain one
        raise : ObservationParsingException
        """
        el = self._get_child_element("visibility", parent, ns, False)
        if el is None:
            return None
        distance = self._get_interval("distance", el, ns, True)
        de = self._get_child_text_as_float("distributionEccentricity", el, ns, True)
        df = self._get_child_text_as_float("distributionFill", el, ns, True)
        return plane.Visibility(distance, de, df)

    def _get_polarization(self, element_tag, parent, ns, required):
        """Build a Polarization object from an XML representation of a
        polarization element.

        Arguments:
        elTag : element tag which identifies the element
        parent : element containing the position element
        ns : namespace of the document
        required : boolean indicating whether the element is required
        return : a Polarization object or
                 None if the document does not contain a polarization element
        raise : ObservationParsingException
        """
        el = self._get_child_element(element_tag, parent, ns, required)
        if el is None:
            return None
        _pstates_el = self._get_child_element("states", el, ns, False)
        if _pstates_el is not None:
            _polarization_states = list()
            for _pstate_el in self._get_children(_pstates_el, ns, "state"):
                if not isinstance(_pstate_el, str):
                    _pstate = _pstate_el.text
                else:
                    _pstate = _pstate_el
                _polarization_state = plane.PolarizationState(_pstate)
                _polarization_states.append(_polarization_state)

            _dimension = self._get_child_text_as_int("dimension", el,
                                                     ns, False)
            return Polarization(dimension=_dimension, states=_polarization_states)
        else:
            return None

    def _get_shape(self, element_tag, parent, ns, required):
        shape_element = self._get_child_element(element_tag, parent, ns,
                                                required)
        if shape_element is None:
            return None
        if self._input_format == 'xml':
            shape_type = shape_element.get(XSI + "type")
        else:
            shape_type = shape_element.get("@type")
        if "caom2:Polygon" == shape_type:
            polygon = self._get_polygon(ns, shape_element)
            samples = None
            if self.version < 25:
                samples_element = self._get_child_element("samples", shape_element,
                                                          ns, True)
                vertices = list()
                self._add_vertices(vertices, samples_element, ns)
                samples = _to_samples(vertices)
            else:
                samples = self._get_child_element("samples", parent,
                                                  ns, True)
                sample_list = []
                for _shape in self._get_children(samples, ns, "shape"):
                    if self._input_format == 'xml':
                        sample_type = _shape.get(XSI + "type")
                    else:
                        sample_type = _shape.get("@type")
                    if "caom2:Polygon" == sample_type:
                        sample_list.append(self._get_polygon(ns, _shape))
                    elif "caom2:Circle" == sample_type:
                        sample_list.append(self._get_circle(ns, _shape))
                    else:
                        raise TypeError("Unsupported sample type " + sample_type)
                samples = shape.MultiShape(sample_list)
            return polygon, samples
        elif "caom2:Circle" == shape_type:
            circle = self._get_circle(ns, shape_element)
            return circle, shape.MultiShape([circle])
        else:
            raise TypeError("Unsupported shape type " + shape_type)

    def _get_circle(self, ns, shape_element):
        center = self._get_child_element("center", shape_element, ns, True)
        center_point = self._get_point(center, ns, True)
        radius = self._get_child_text_as_float(
            "radius", shape_element, ns, True)
        circle = shape.Circle(center=center_point, radius=radius)
        return circle

    def _get_polygon(self, ns, shape_element):
        points_element = self._get_child_element("points", shape_element,
                                                 ns, True)
        points = list()
        for point in self._get_children(points_element, ns, "point"):
            points.append(self._get_point(point, ns, True))
        return shape.Polygon(points)

    def _add_energy_bands(self, energy_bands, parent, ns):
        """Create EnergyBand objects from an XML representation of
        observation URI elements found in energy_band element, and add them to
        the set of energy_bads

        Arguments:
        energy_bands : Set of energy bands to add to
        parent : element containing the Environment element
        ns : namespace of the document
        raise : ObservationParsingException
        """
        # this is for backwards compatibility with pre 2.4 versions
        em_band = self._get_child_text("emBand", parent, ns, False)
        if em_band:
            energy_bands.add(plane.EnergyBand(em_band))

        el = self._get_child_element("energyBands", parent, ns, False)
        if el is not None:
            for eb in self._get_children(el, ns, "emBand"):
                if not isinstance(eb, str):
                    eb = eb.text
                energy_bands.add(plane.EnergyBand(eb))

    def _add_vertices(self, vertices, parent, ns):
        _vertices_element = self._get_child_element("vertices", parent, ns,
                                                    False)
        if _vertices_element is None:
            return None
        else:
            for _vertex_element in self._get_children(_vertices_element, ns, "vertex"):
                cval1 = self._get_child_text_as_float("cval1", _vertex_element,
                                                      ns, False)
                cval2 = self._get_child_text_as_float("cval2", _vertex_element,
                                                      ns, False)
                seg_type_value = self._get_child_text_as_int("type",
                                                             _vertex_element,
                                                             ns, False)
                seg_type = shape.SegmentType(seg_type_value)
                _vertex = shape.Vertex(cval1, cval2, seg_type)
                vertices.append(_vertex)

    def _get_interval(self, element_tag, parent, ns, required):
        _interval_el = self._get_child_element(element_tag, parent, ns,
                                               required)
        if _interval_el is None:
            return None
        _lower = self._get_child_text_as_float("lower", _interval_el, ns, True)
        _upper = self._get_child_text_as_float("upper", _interval_el, ns, True)
        _interval = shape.Interval(_lower, _upper)
        return _interval

    def _get_samples(self, parent, ns, required):
        _samples_el = self._get_child_element("samples", parent, ns, required)
        if _samples_el is None:
            return None
        _samples = []
        for _sample_el in self._get_children(_samples_el, ns, "sample"):
            _si_lower = self._get_child_text_as_float("lower", _sample_el,
                                                      ns, required)
            _si_upper = self._get_child_text_as_float("upper", _sample_el,
                                                      ns, required)
            _sub_interval = dali.Interval(_si_lower, _si_upper)
            _samples.append(_sub_interval)
        return _samples

    def _add_chunks(self, chunks, parent, ns):
        """Build Chunk objects from an XML representation of Chunk elements
        and add them to the set of Chunks.

        Argument:
        chunks : set of Chunk objects from the Part
        parent : element containing the Chunk elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._get_child_element("chunks", parent, ns, False)
        if el is None:
            return None
        else:
            for chunk_element in self._get_children(el, ns, "chunk"):
                _chunk = chunk.Chunk()
                product_type = \
                    self._get_child_text("productType", chunk_element, ns,
                                         False)
                if product_type:
                    _chunk.product_type = \
                        chunk.DataLinkSemantics(product_type)
                _chunk.naxis = \
                    self._get_child_text_as_int("naxis", chunk_element, ns,
                                                False)
                _chunk.observable_axis = \
                    self._get_child_text_as_int("observableAxis",
                                                chunk_element, ns, False)
                _chunk.position_axis_1 = \
                    self._get_child_text_as_int("positionAxis1", chunk_element,
                                                ns, False)
                _chunk.position_axis_2 = \
                    self._get_child_text_as_int("positionAxis2", chunk_element,
                                                ns, False)
                _chunk.energy_axis = \
                    self._get_child_text_as_int("energyAxis", chunk_element,
                                                ns, False)
                _chunk.time_axis = \
                    self._get_child_text_as_int("timeAxis", chunk_element, ns,
                                                False)
                _chunk.custom_axis = \
                    self._get_child_text_as_int("customAxis", chunk_element,
                                                ns, False)
                _chunk.polarization_axis = \
                    self._get_child_text_as_int("polarizationAxis",
                                                chunk_element, ns, False)
                _chunk.observable = \
                    self._get_observable_axis("observable", chunk_element, ns,
                                              False)
                _chunk.position = \
                    self._get_spatial_wcs("position", chunk_element, ns, False)
                _chunk.energy = \
                    self._get_spectral_wcs("energy", chunk_element, ns, False)
                _chunk.time = \
                    self._get_temporal_wcs("time", chunk_element, ns, False)
                _chunk.polarization = \
                    self._get_polarization_wcs("polarization", chunk_element,
                                               ns, False)
                _chunk.custom = \
                    self._get_custom_wcs("custom", chunk_element,
                                         ns, False)
                self._set_entity_attributes(chunk_element, ns, _chunk)
                chunks.append(_chunk)

    def _add_parts(self, parts, parent, ns):
        """Build Part objects from an XML representation of Part elements and
        add them to the set of Parts.

        Argument:
        parts : set of Part objects from the Artifact
        parent : element containing the Part elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._get_child_element("parts", parent, ns, False)
        if el is None:
            return None
        else:
            for part_element in self._get_children(el, ns, "part"):
                _part = \
                    part.Part(
                        self._get_child_text("name", part_element, ns, True))
                product_type = \
                    self._get_child_text("productType", part_element, ns,
                                         False)
                if product_type:
                    _part.product_type = \
                        chunk.DataLinkSemantics(product_type)
                self._add_chunks(_part.chunks, part_element, ns)
                self._set_entity_attributes(part_element, ns, _part)
                parts[_part.name] = _part

    def _add_artifacts(self, artifacts, parent, ns):
        """Build artifacts from an XML representation of the artifact elements
        and add them to the set of Artifacts.

        Arguments:
        artifacts : set of Artifacts from the Plane
        parent : element containing the Artifact elements
        ns : namespace fo the document
        raise : ObservationParsingException
        """
        el = self._get_child_element("artifacts", parent, ns, False)
        if el is None:
            return None
        else:
            for artifact_element in self._get_children(el, ns, "artifact"):
                uri = self._get_child_text("uri", artifact_element, ns, True)

                product_type = self._get_child_text("productType",
                                                    artifact_element, ns,
                                                    False)
                if product_type is None:
                    product_type = chunk.DataLinkSemantics.SCIENCE
                    print(
                        "Using default Artifact.productType value {0}".format(
                            str(chunk.DataLinkSemantics.SCIENCE)))
                else:
                    product_type = chunk.DataLinkSemantics(product_type)

                release_type = self._get_child_text("releaseType",
                                                    artifact_element, ns,
                                                    False)
                if release_type is None:
                    release_type = artifact.ReleaseType.DATA
                    print(
                        "Using default Artifact.releaseType value {0}".format(
                            str(artifact.ReleaseType.DATA)))
                else:
                    release_type = artifact.ReleaseType(release_type)

                _artifact = artifact.Artifact(uri, product_type, release_type)
                if self.version >= 25:
                    sub = self._get_child_text("uriBucket", artifact_element, ns, True)
                    if sub != _artifact.uri_bucket:
                        raise ObservationParsingException(
                            "Parsed artifact URI bucket {} does not match calculated artifact URI bucket {}".
                            format(sub, _artifact.uri_bucket))
                _artifact.description_id = self._get_child_text("descriptionID",
                                                                artifact_element,
                                                                ns, False)
                cr = self._get_child_text("contentRelease", artifact_element,
                                          ns, False)
                _artifact.content_release = caom_util.str2ivoa(cr)
                _artifact.content_read_groups = \
                    self._get_groups("contentReadGroups", artifact_element,
                                     ns, False)
                _artifact.content_type = self._get_child_text("contentType",
                                                              artifact_element,
                                                              ns, False)
                _artifact.content_length = (
                    self._get_child_text_as_long("contentLength",
                                                 artifact_element, ns, False))
                content_checksum = self._get_child_text("contentChecksum",
                                                        artifact_element, ns,
                                                        False)
                if content_checksum:
                    _artifact.content_checksum = content_checksum
                self._add_parts(_artifact.parts, artifact_element, ns)
                self._set_entity_attributes(artifact_element, ns, _artifact)
                artifacts[_artifact.uri] = _artifact

    def _get_children(self, parent, ns, child_tag):
        """Get a list of child elements. For XML input, children have a tag.
        For JSON input, children are in a list in the parent."""
        if self._input_format == 'xml':
            plane_list = parent.iterchildren('{' + ns + '}' + child_tag)
        else:
            plane_list = parent
        return plane_list

    def _add_planes(self, obs, parent, ns):
        """Create Planes object from XML representation of Plane elements
        and add them to the set of Planes.

        Arguments:
        obs : Observation object containing the Planes
        parent : element containing the Plane elements
        ns : namespace of the document
        raise : ObservationParsingException
        """
        el = self._get_child_element("planes", parent, ns, False)
        if el is None:
            return None
        else:
            for plane_element in self._get_children(el, ns, "plane"):
                if self.version < 25:
                    # TODO
                    _uri = "{}/{}".format(
                        obs.uri,
                        self._get_child_text("productID",
                                             plane_element, ns, True))
                else:
                    _uri = self._get_child_text("uri", plane_element, ns, True)
                _plane = plane.Plane(_uri)
                _plane.meta_release = caom_util.str2ivoa(
                    self._get_child_text("metaRelease", plane_element, ns,
                                         False))
                _plane.data_release = caom_util.str2ivoa(
                    self._get_child_text("dataRelease", plane_element, ns,
                                         False))
                _plane.meta_read_groups = self._get_groups("metaReadGroups",
                                                           plane_element, ns,
                                                           False)
                _plane.data_read_groups = self._get_groups("dataReadGroups",
                                                           plane_element, ns,
                                                           False)
                data_product_type = \
                    self._get_child_text("dataProductType", plane_element, ns,
                                         False)
                if data_product_type:
                    _plane.data_product_type = plane.DataProductType(data_product_type)
                calibration_level = \
                    self._get_child_text("calibrationLevel", plane_element, ns,
                                         False)
                if calibration_level:
                    _plane.calibration_level = \
                        plane.CalibrationLevel(int(calibration_level))
                _plane.provenance = \
                    self._get_provenance("provenance", plane_element, ns,
                                         False)
                _plane.observable = self._get_observable(plane_element, ns)

                _plane.metrics = \
                    self._get_metrics("metrics", plane_element, ns, False)
                _plane.quality = \
                    self._get_quality("quality", plane_element, ns, False)

                _plane.position = self._get_position("position", plane_element,
                                                     ns, False)
                _plane.energy = self._get_energy("energy", plane_element, ns,
                                                 False)
                _plane.time = self._get_time("time", plane_element, ns, False)
                _plane.polarization = self._get_polarization("polarization",
                                                             plane_element, ns,
                                                             False)
                _plane.custom = \
                    self._get_custom("custom", plane_element, ns, False)
                _plane.visibility = self._get_visibility(plane_element, ns)
                self._add_artifacts(_plane.artifacts, plane_element, ns)
                self._set_entity_attributes(plane_element, ns, _plane)
                obs.planes[_plane.uri] = _plane

    def read(self, source):
        """Build an Observation object from an XML document located in source.
        Source an be a file name/path, a file object, a file-like object or a
        URL using the HTTP or FTP protocol.

        Arguments:
        source : source of XML document containing an Observation element
        return : an Observation object
        raise : ObservationParsingException
        """

        if self._input_format == 'xml':
            doc = etree.parse(source)
            if self._validate and self._xmlschema:
                self._xmlschema.assertValid(doc)
            root = doc.getroot()
            ns = root.nsmap["caom2"]
            self.version = CAOM_VERSION[ns]
        else:  # json
            src = source
            if isinstance(source, IOBase):
                source.seek(0)  # Ensure we're reading from the start
                src = source.read()

            # Case 2: source is a path to an existing file
            elif isinstance(source, str) and os.path.isfile(source):
                with open(source, 'r', encoding='utf-8') as f:
                    src = f.read()

            if not isinstance(src, str):
                raise TypeError(
                    "Unsupported source type. Must be string, stream, or file path.")
            root = json.loads(src)
            if "http://www.opencadc.org/caom2/xml/v2.5" == root.get("@caom2", None):
                self.version = 25
            else:
                self.version = 24
            ns = None
        collection = self._get_child_text("collection", root, ns, True)
        if self._input_format == 'xml' and self.version < 25:
            observation_id = \
                str(self._get_child_text("observationID", root, ns, True))
            uri = "caom:" + collection + "/" + observation_id
        else:
            uri = self._get_child_text("uri", root, ns, True)

        # Instantiate Algorithm
        algorithm = self._get_algorithm("algorithm", root, ns, True)
        # Instantiate Observation
        if root.get("{http://www.w3.org/2001/XMLSchema-instance}type") \
                == "caom2:SimpleObservation":
            obs = observation.SimpleObservation(collection, uri)
            obs.algorithm = algorithm
        else:
            obs = \
                observation.DerivedObservation(collection, uri, algorithm)
        if self.version >= 25:
            sub = self._get_child_text("uriBucket", root, ns, True)
            if sub != obs.uri_bucket:
                raise ObservationParsingException(
                    "Parsed obs URI bucket {} does not match calculated obs URI bucket {}".
                    format(sub, obs.uri_bucket))
        # Instantiate children of Observation
        obs.sequence_number = \
            self._get_child_text_as_int("sequenceNumber", root, ns, False)
        obs.type = \
            self._get_child_text("type", root, ns, False)
        intent = self._get_child_text("intent", root, ns, False)
        if intent:
            obs.intent = observation.ObservationIntentType(intent)
        obs.meta_release = \
            self._get_meta_release("metaRelease", root, ns, False)
        obs.meta_read_groups = \
            self._get_groups('metaReadGroups', root, ns, False)
        obs.proposal = \
            self._get_proposal("proposal", root, ns, False)
        obs.target = \
            self._get_target("target", root, ns, False)
        obs.target_position = \
            self._get_target_position("targetPosition", root, ns, False)
        obs.telescope = \
            self._get_telescope("telescope", root, ns, False)
        obs.instrument = \
            self._get_instrument("instrument", root, ns, False)
        obs.environment = \
            self._get_environment("environment", root, ns, False)
        obs.requirements = \
            self._get_requirements("requirements", root, ns, False)
        self._add_planes(obs, root, ns)
        if isinstance(obs, observation.DerivedObservation):
            self._add_members(obs.members, root, ns)

        self._set_entity_attributes(root, ns, obs)
        return obs


class ObservationWriter(object):
    """ ObservationWriter """

    def __init__(self, validate=False, write_empty_collections=False,
                 namespace_prefix="caom2", namespace=None, output_format='xml'):
        """
        Arguments:
        validate : If True enable schema validation, False otherwise
        write_empty_collections : if True write empty elements for empty
        collections namespace_prefix : a CAOM-2.x namespace prefix
        namespace : a valid CAOM-2.x target namespace
        """
        if output_format not in ['xml', 'json']:
            raise RuntimeError('invalid output format {}'.format(output_format))
        self._output_format = output_format

        self._validate = validate
        self._write_empty_collections = write_empty_collections

        if (namespace_prefix is None or not namespace_prefix) \
                and output_format == 'xml':
            raise RuntimeError('null or empty namespace_prefix not allowed')

        if namespace is None or namespace == CAOM25_NAMESPACE:
            self._output_version = 25
            self._caom2_namespace = CAOM25
            self._namespace = CAOM25_NAMESPACE
        elif namespace == CAOM24_NAMESPACE:
            self._output_version = 24
            self._caom2_namespace = CAOM24
            self._namespace = CAOM24_NAMESPACE
        elif namespace == CAOM23_NAMESPACE:
            self._output_version = 23
            self._caom2_namespace = CAOM23
            self._namespace = CAOM23_NAMESPACE
        else:
            raise RuntimeError('invalid namespace {}'.format(namespace))

        if self._validate:
            if self._output_format != 'xml':
                raise AttributeError("Validation works with xml output only")
            if self._output_version == 25:
                schema_file = CAOM25_SCHEMA_FILE
            elif self._output_version == 24:
                schema_file = CAOM24_SCHEMA_FILE
            elif self._output_version == 23:
                schema_file = CAOM23_SCHEMA_FILE
            schema_path = os.path.join(THIS_DIR + '/' + DATA_PKG,
                                       schema_file)
            # schema_path = pkg_resources.resource_filename(
            #     DATA_PKG, schema_file)
            xmlschema_doc = etree.parse(schema_path)
            self._xmlschema = etree.XMLSchema(xmlschema_doc)

        self._nsmap = {namespace_prefix: self._namespace, "xsi": XSI_NAMESPACE}

    def _to_json(tree):
        def xml_to_dict(element):
            node = {}

            # Add attributes if any
            for attrib in element.attrib:
                localname = etree.QName(attrib).localname
                node['@' + localname] = element.attrib[attrib]

            # Add children
            children = list(element)
            if children:
                child_dict = {}
                for child in children:
                    child_tag = etree.QName(child).localname
                    child_data = xml_to_dict(child)
                    if child_tag in child_dict:
                        # Convert to list if multiple children with same tag
                        if not isinstance(child_dict[child_tag], list):
                            child_dict[child_tag] = [child_dict[child_tag]]
                        child_dict[child_tag].append(child_data)
                    else:
                        child_dict[child_tag] = child_data
                keys = list(child_dict.keys())
                if len(keys) == 1 and isinstance(child_dict[keys[0]], list):
                    node = child_dict[keys[0]]
                else:
                    node.update(child_dict)
            else:
                # No children — use text
                text = element.text.strip() if element.text else ''
                node = text if not element.attrib else {**node, '#text': text}

            return node

        # Convert to dict
        data_dict = xml_to_dict(tree)

        # Convert dict to JSON
        return json.dumps(data_dict, indent=2)

    def write(self, obs, out):
        """
        Writes an observation to an output stream
        :param obs: observation to write
        :param out: stream or file name to write to
        :return:
        """
        caom_util.type_check(obs, observation.Observation, "observation")

        if isinstance(obs, observation.SimpleObservation):
            obs_type = "caom2:SimpleObservation"
        else:
            if self._output_version < 24:
                obs_type = "caom2:CompositeObservation"
            else:
                obs_type = "caom2:DerivedObservation"
        if self._output_format == 'xml':
            obs_element = etree.Element(self._caom2_namespace + "Observation",
                                        nsmap=self._nsmap)
            obs_element.set(XSI + "type", obs_type)
        else:
            if self._output_version == 25:
                caom2_v = 'v2.5'
            else:
                caom2_v = 'v2.4'
            obs_element = {'@caom2': "http://www.opencadc.org/caom2/xml/" + caom2_v,
                           '@type': obs_type}

        self._add_entity_attributes(obs, obs_element)

        self._add_element("collection", obs.collection, obs_element)
        if self._output_version < 25:
            observation_id = obs.uri.split('/')[-1]
            self._add_element("observationID", observation_id, obs_element)
        else:
            self._add_element("uri", obs.uri, obs_element)
            self._add_element('uriBucket', obs.uri_bucket, obs_element)
        self._add_datetime_element("metaRelease", obs.meta_release,
                                   obs_element)
        if self._output_version < 24 and obs.meta_read_groups:
            raise AttributeError(
                "Attempt to write CAOM2.4 element "
                "(observation.meta_read_groups) as CAOM2.3 Observation")
        self._add_groups_element("metaReadGroups", obs.meta_read_groups,
                                 obs_element)
        self._add_element("sequenceNumber", obs.sequence_number, obs_element)
        self._add_algorithm_element(obs.algorithm, obs_element)
        self._add_element("type", obs.type, obs_element)
        if obs.intent is not None:
            self._add_element(
                "intent", obs.intent.value, obs_element)

        self._add_proposal_element(obs.proposal, obs_element)
        self._add_target_element(obs.target, obs_element)
        self._add_target_position_element(obs.target_position, obs_element)
        self._add_requirements_element(obs.requirements, obs_element)
        self._add_telescope_element(obs.telescope, obs_element)
        self._add_instrument_element(obs.instrument, obs_element)
        self._add_environment_element(obs.environment, obs_element)
        self._add_planes_element(obs.planes, obs_element)

        if isinstance(obs, observation.DerivedObservation):
            self._add_members_element(obs.members, obs_element)

        if self._validate and self._xmlschema:
            self._xmlschema.assertValid(obs_element)

        if self._output_format == 'xml':
            tree = etree.ElementTree(obs_element)
            try:
                # try to write as text first
                out.write(etree.tostring(tree, encoding='unicode',
                                         pretty_print=True))
                return
            except Exception:
                pass  # didn't work. Try to write as binary
            tree.write(out, encoding='utf-8',
                       xml_declaration=True, pretty_print=True)
        else:
            out.write(json.dumps(obs_element, indent=2, separators=(',', ' : ')))

    def _add_entity_attributes(self, entity, element):
        self._add_attribute("id", str(entity._id), element)

        if entity._last_modified is not None:
            self._add_attribute(
                "lastModified", caom_util.date2ivoa(entity._last_modified),
                element)

        if entity._max_last_modified is not None:
            self._add_attribute(
                "maxLastModified",
                caom_util.date2ivoa(entity._max_last_modified), element)
        if entity._meta_checksum is not None:
            self._add_attribute(
                "metaChecksum", entity._meta_checksum, element)
        if entity._acc_meta_checksum is not None:
            self._add_attribute(
                "accMetaChecksum", entity._acc_meta_checksum, element)

        if self._output_version >= 24:
            if entity._meta_producer is not None:
                self._add_attribute(
                    "metaProducer", entity.meta_producer, element)

    def _add_algorithm_element(self, algorithm, parent):
        if algorithm is None:
            return

        element = self._get_caom_element("algorithm", parent)
        self._add_element("name", algorithm.name, element)

    def _add_proposal_element(self, proposal, parent):
        if proposal is None:
            return

        element = self._get_caom_element("proposal", parent)
        self._add_element("id", proposal.id, element)
        self._add_element("pi", proposal.pi, element)
        self._add_element("project", proposal.project, element)
        self._add_element("title", proposal.title, element)
        self._add_element("reference", proposal.reference, element)
        self._add_keywords_element(proposal.keywords, element)

    def _add_target_element(self, target, parent):
        if target is None:
            return

        element = self._get_caom_element("target", parent)
        self._add_element("name", target.name, element)
        if target.target_id is not None:
            if self._output_version >= 24:
                self._add_element("targetID", target.target_id, element)
            else:
                raise AttributeError(
                    "Attempt to write CAOM2.4 element (target.targetID) "
                    "as CAOM2.3 Observation")
        if target.type is not None:
            self._add_element("type", target.type.value, element)
        self._add_boolean_element("standard", target.standard, element)
        self._add_element("redshift", target.redshift, element)
        self._add_boolean_element("moving", target.moving, element)
        self._add_keywords_element(target.keywords, element)

    def _add_target_position_element(self, target_position, parent):
        if target_position is None:
            return

        element = self._get_caom_element("targetPosition", parent)
        self._add_element("coordsys", target_position.coordsys, element)
        if target_position.equinox is not None:
            self._add_element("equinox", target_position.equinox, element)
        self._add_point_element("coordinates", target_position.coordinates,
                                element)

    def _add_requirements_element(self, requirements, parent):
        if requirements is None:
            return

        element = self._get_caom_element("requirements", parent)
        self._add_element(
            "flag", requirements.flag.value, element)

    def _add_telescope_element(self, telescope, parent):
        if telescope is None:
            return

        element = self._get_caom_element("telescope", parent)
        self._add_element("name", telescope.name, element)
        self._add_element("geoLocationX", telescope.geo_location_x, element)
        self._add_element("geoLocationY", telescope.geo_location_y, element)
        self._add_element("geoLocationZ", telescope.geo_location_z, element)
        if telescope.tracking_mode:
            self._add_element("trackingMode", telescope.tracking_mode.value, element)
        self._add_keywords_element(telescope.keywords, element)

    def _add_instrument_element(self, instrument, parent):
        if instrument is None:
            return

        element = self._get_caom_element("instrument", parent)
        self._add_element("name", instrument.name, element)
        self._add_keywords_element(instrument.keywords, element)

    def _add_environment_element(self, environment, parent):
        if environment is None:
            return

        element = self._get_caom_element("environment", parent)
        self._add_element("seeing", environment.seeing, element)
        self._add_element("humidity", environment.humidity, element)
        self._add_element("elevation", environment.elevation, element)
        self._add_element("tau", environment.tau, element)
        self._add_element("wavelengthTau", environment.wavelength_tau, element)
        self._add_element("ambientTemp", environment.ambient_temp, element)
        self._add_boolean_element("photometric", environment.photometric,
                                  element)

    def _add_members_element(self, members, parent):
        if members is None or \
                (len(members) == 0 and not self._write_empty_collections):
            return

        element = self._get_caom_element("members", parent)
        if self._output_format == 'xml':
            if self._output_version < 25:
                member_tag = "observationURI"
            else:
                member_tag = "member"
            for member in members:
                member_element = self._get_caom_element(member_tag, element)
                member_element.text = member
        else:
            parent["members"] = list(members) if members else []

    def _add_groups_element(self, name, groups, parent):
        if self._output_version < 24:
            return
        if not groups and not self._write_empty_collections:
            return
        element = self._get_caom_element(name, parent)
        if self._output_format == 'xml':
            for gr in groups:
                gr_elem = self._get_caom_element("groupURI", element)
                gr_elem.text = gr
        else:
            parent[name] = sorted(list(groups)) if groups else []

    def _add_planes_element(self, planes, parent):
        if planes is None or \
                (len(planes) == 0 and not self._write_empty_collections):
            return
        if self._output_format == 'xml':
            element = self._get_caom_element("planes", parent)
        else:
            element = parent["planes"] = []
        for _plane in planes.values():
            if self._output_format == 'xml':
                plane_element = self._get_caom_element("plane", element)
            else:
                plane_element = {}
                element.append(plane_element)
            self._add_entity_attributes(_plane, plane_element)
            if self._output_version < 25:
                _comp = _plane.uri.split('/')
                if len(_comp) != 3:
                    raise ValueError("Attempt to write CAOM2.4 but can't deduce "
                                     "Plane.productID in Plane.uri=" + _plane.uri)
                self._add_element("productID", _comp[-1], plane_element)
            else:
                self._add_element("uri", _plane.uri, plane_element)
            self._add_datetime_element("metaRelease", _plane.meta_release,
                                       plane_element)
            if self._output_version < 24 and _plane.meta_read_groups:
                raise AttributeError(
                    "Attempt to write CAOM2.4 element "
                    "(plane.meta_read_groups) as CAOM2.3 Observation")
            self._add_groups_element("metaReadGroups", _plane.meta_read_groups,
                                     plane_element)
            self._add_datetime_element("dataRelease", _plane.data_release,
                                       plane_element)
            if self._output_version < 24 and _plane.data_read_groups:
                raise AttributeError(
                    "Attempt to write CAOM2.4 element "
                    "(plane.data_read_groups) as CAOM2.3 Observation")
            self._add_groups_element("dataReadGroups", _plane.data_read_groups,
                                     plane_element)
            if _plane.data_product_type is not None:
                self._add_element("dataProductType",
                                  _plane.data_product_type.value,
                                  plane_element)
            if _plane.calibration_level is not None:
                self._add_element("calibrationLevel",
                                  _plane.calibration_level.value,
                                  plane_element)
            self._add_provenance_element(_plane.provenance, plane_element)
            self._add_observable_element(_plane.observable, plane_element)
            self._add_metrics_element(_plane.metrics, plane_element)
            self._add_quality_element(_plane.quality, plane_element)

            self._add_position_element(_plane.position, plane_element)
            self._add_energy_element(_plane.energy, plane_element)
            self._add_time_element(_plane.time, plane_element)
            self._add_polarization_element(_plane.polarization,
                                           plane_element)
            if _plane.custom is not None:
                if self._output_version < 24:
                    raise AttributeError(
                        'plane custom element only available in CAOM2.4')
                else:
                    self._add_custom_element(_plane.custom, plane_element)

            self._add_visibility_element(_plane.visibility, plane_element)

            self._add_artifacts_element(_plane.artifacts, plane_element)

    def _add_visibility_element(self, visibility, parent):
        if visibility is None:
            return

        if self._output_version < 25:
            raise AttributeError("Attempt to output CAOM2.5 attribute (Plane.visibility) as "
                                 "{} Observation".format(self._output_version))

        element = self._get_caom_element("visibility", parent)
        self._add_interval_element("distance", visibility.distance, element)
        self._add_element("distributionEccentricity", visibility.distribution_eccentricity,
                          element)
        self._add_element("distributionFill", visibility.distribution_fill, element)

    def _add_position_element(self, position, parent):
        if position is None:
            return

        element = self._get_caom_element("position", parent)
        self._add_bounds_and_samples(position, element)
        if position.min_bounds:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (position.min_bounds) as "
                    "{} Observation".format(self._output_version))
            else:
                self._add_shape_element("minBounds", element, position.min_bounds)
        self._add_dimension2d_element("dimension", position.dimension, element)
        if position.max_recoverable_scale:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (position.max_recoverable_scale) as "
                    "{} Observation".format(self._output_version))
            else:
                self._add_interval_element("maxRecoverableScale", position.max_recoverable_scale, element)
        self._add_element("resolution", position.resolution, element)
        if self._output_version < 24:
            if position.resolution_bounds is not None:
                raise AttributeError(
                    "Attempt to write CAOM2.4 element "
                    "(position.resolution_bounds) as CAOM2.3 Observation")
        else:
            self._add_interval_element("resolutionBounds",
                                       position.resolution_bounds, element)
        self._add_element("sampleSize", position.sample_size, element)
        if position.calibration:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (position.calibration) as "
                    "{} Observation".format(self._output_version))
            else:
                self._add_element("calibration", position.calibration.value, element)

    def _add_energy_element(self, energy, parent):
        if energy is None:
            return
        element = self._get_caom_element("energy", parent)
        bounds_elem = self._add_interval_element("bounds", energy.bounds, element)
        if self._output_version < 25:
            # samples element is within bounds element
            self._add_samples_element(energy.samples, bounds_elem)
        else:
            self._add_samples_element(energy.samples, element)

        self._add_element("dimension", energy.dimension, element)
        self._add_element("resolvingPower", energy.resolving_power, element)
        if energy.resolving_power_bounds is not None:
            if self._output_version < 24:
                if len(energy.energy_bands) > 1:
                    raise AttributeError(
                        "Attempt to write CAOM2.4 element "
                        "(energy.resolving_power_bands) as "
                        "CAOM2.3 Observation")
            else:
                self._add_interval_element("resolvingPowerBounds",
                                           energy.resolving_power_bounds, element)
        if energy.resolution is not None:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (energy.resolution) as "
                    "{} Observation".format(self._output_version))
            else:
                self._add_element("resolution", energy.resolution, element)
        if energy.resolution_bounds is not None:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (energy.resolution_bounds) as "
                    "{} Observation".format(self._output_version))
            else:
                self._add_interval_element("resolutionBounds",
                                           energy.resolution_bounds, element)

        self._add_element("sampleSize", energy.sample_size, element)
        self._add_element("bandpassName", energy.bandpass_name, element)
        if energy.energy_bands:
            if self._output_version < 24:
                if len(energy.energy_bands) > 1:
                    raise AttributeError(
                        "Attempt to write CAOM2.4 element "
                        "(energy.energy_bands) as CAOM2.3 Observation")
                else:
                    self._add_element("emBand",
                                      next(iter(energy.energy_bands)).value,
                                      element)
            else:
                self._add_element("energyBands", None, element)
                if self._output_format == 'xml':
                    eb_element = self._get_caom_element("energyBands", element)
                else:
                    eb_element = element["energyBands"] = []
                for bb in sorted(energy.energy_bands):
                    self._add_element("emBand", bb.value, eb_element)
        if self._output_version < 25:
            self._add_element("restwav", energy.rest, element)
        else:
            self._add_element("rest", energy.rest, element)
        if energy.transition:
            transition = self._get_caom_element("transition", element)
            self._add_element("species", energy.transition.species, transition)
            self._add_element("transition", energy.transition.transition,
                              transition)
        if energy.calibration:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (energy.calibration) as "
                    "{} Observation".format(self._output_version))
            else:
                self._add_element("calibration", energy.calibration.value, element)

    def _add_time_element(self, time, parent):
        if time is None:
            return
        element = self._get_caom_element("time", parent)
        bounds_elem = self._add_interval_element("bounds", time.bounds, element)
        if self._output_version < 25:
            # samples element is within bounds element
            self._add_samples_element(time.samples, bounds_elem)
        else:
            self._add_samples_element(time.samples, element)
        self._add_element("dimension", time.dimension, element)
        self._add_element("resolution", time.resolution, element)
        if self._output_version < 24:
            if time.resolution_bounds is not None:
                raise AttributeError(
                    "Attempt to write CAOM2.4 element "
                    "(time.resolution_bounds) as CAOM2.3 Observation")
        else:
            self._add_interval_element("resolutionBounds",
                                       time.resolution_bounds, element)
        self._add_element("sampleSize", time.sample_size, element)
        self._add_element("exposure", time.exposure, element)
        if time.exposure_bounds is not None:
            if self._output_version < 25:
                if len(time.exposure_bounds) > 1:
                    raise AttributeError(
                        "Attempt to write CAOM2.5 element "
                        "(time.exposure_bounds) as {} Observation".format(self._output_version))
            else:
                self._add_interval_element("exposureBounds", time.exposure_bounds, element)
        if time.calibration:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (time.calibration) as {} "
                    "Observation".format(self._output_version))
            else:
                self._add_element("calibration", time.calibration.value, element)

    def _add_custom_element(self, custom, parent):
        if custom is None:
            return
        element = self._get_caom_element("custom", parent)
        self._add_element("ctype", custom.ctype, element)
        bounds_elem = self._add_interval_element("bounds", custom.bounds, element)
        if self._output_version < 25:
            # samples element is within bounds element
            self._add_samples_element(custom.samples, bounds_elem)
        else:
            self._add_samples_element(custom.samples, element)
        self._add_element("dimension", custom.dimension, element)

    def _add_polarization_element(self, polarization, parent):
        if polarization is None:
            return
        element = self._get_caom_element("polarization", parent)
        if not polarization.states:
            raise AttributeError("Polarization.states missing")
        if polarization.states:
            if self._output_format == 'xml':
                _pstates_el = self._get_caom_element("states", element)
            else:
                _pstates_el = element["states"] = []
            for _state in polarization.states:
                self._add_element("state", _state.value, _pstates_el)
        self._add_element("dimension", polarization.dimension, element)

    def _add_bounds_and_samples(self, position, parent):
        if position is None:
            return

        bounds_element = self._add_shape_element("bounds", parent, position.bounds)
        if self._output_version < 25:
            if isinstance(position.bounds, shape.Circle):
                # samples not supported for circles so need to make sure that
                # samples is the same with bounds
                if len(position.samples.shapes) > 1 or position.samples.shapes[0] != position.bounds:
                    raise AttributeError(
                        "Cannot write a CAOM25 circle bounds position as CAOM{} "
                        "Observation".format(self._output_version))
                else:
                    pass
            else:
                samples_element = self._get_caom_element("samples", bounds_element)
                vertices_element = self._get_caom_element("vertices",
                                                          samples_element)
                vertices = _to_vertices(position.samples)
                for vertex in vertices:
                    vertex_element = self._get_caom_element("vertex",
                                                            vertices_element)
                    self._add_element("cval1", vertex.cval1, vertex_element)
                    self._add_element("cval2", vertex.cval2, vertex_element)
                    self._add_element("type", vertex.type.value, vertex_element)
        else:
            if self._output_format == 'xml':
                samples_element = self._get_caom_element("samples", parent=parent)
            else:
                samples_element = parent["samples"] = []
            for samples_shape in position.samples.shapes:
                self._add_shape_element("shape", samples_element, samples_shape)

    def _add_shape_element(self, name, parent, elem_shape):
        if isinstance(elem_shape, shape.Polygon):
            if self._output_format == 'xml':
                shape_element = self._get_caom_element(name, parent)
                shape_element.set(XSI + "type", "caom2:Polygon")
                points_element = self._get_caom_element("points",
                                                        shape_element)
            else:
                shape_element = {}
                if isinstance(parent, list):
                    parent.append(shape_element)
                else:
                    parent[name] = shape_element
                shape_element["@type"] = "caom2:Polygon"
                points_element = shape_element["points"] = []
            for point in elem_shape.points:
                self._add_point_element("point", point, points_element)
        elif isinstance(elem_shape, shape.Circle):
            if self._output_format == 'xml':
                shape_element = self._get_caom_element(name, parent)
                shape_element.set(XSI + "type", "caom2:Circle")
            else:
                shape_element = {}
                if isinstance(parent, list):
                    parent.append(shape_element)
                else:
                    parent[name] = shape_element
                shape_element["@type"] = "caom2:Circle"
            self._add_point_element("center", elem_shape.center,
                                    shape_element)
            self._add_element("radius", elem_shape.radius, shape_element)
        else:
            raise TypeError("Unsupported shape type " + elem_shape.__class__.__name__)
        return shape_element

    def _add_interval_element(self, name, interval, parent):
        if interval is None:
            return
        _interval_element = self._get_caom_element(name, parent)
        self._add_element("lower", interval.lower, _interval_element)
        self._add_element("upper", interval.upper, _interval_element)
        return _interval_element

    def _add_samples_element(self, samples, parent):
        if not samples:
            raise AttributeError("non empty samples attribute is required")

        if self._output_format == 'xml':
            _samples_element = self._get_caom_element("samples", parent)
        else:
            _samples_element = parent["samples"] = []
        for _sample in samples:
            if self._output_format == 'xml':
                _sample_element = self._get_caom_element("sample",
                                                         _samples_element)
            else:
                _sample_element = {}
                _samples_element.append(_sample_element)
            self._add_element("lower", _sample.lower, _sample_element)
            self._add_element("upper", _sample.upper, _sample_element)

    def _add_provenance_element(self, provenance, parent):
        if provenance is None:
            return

        element = self._get_caom_element("provenance", parent)
        self._add_element("name", provenance.name, element)
        self._add_element("version", provenance.version, element)
        self._add_element("project", provenance.project, element)
        self._add_element("producer", provenance.producer, element)
        self._add_element("runID", provenance.run_id, element)
        self._add_element("reference", provenance.reference, element)
        self._add_datetime_element("lastExecuted", provenance.last_executed,
                                   element)
        self._add_keywords_element(provenance.keywords, element)
        self._add_inputs_element("inputs", provenance.inputs, element)

    def _add_metrics_element(self, metrics, parent):
        if metrics is None:
            return

        element = self._get_caom_element("metrics", parent)
        self._add_element("sourceNumberDensity", metrics.source_number_density,
                          element)
        self._add_element("background", metrics.background, element)
        self._add_element("backgroundStddev", metrics.background_std_dev,
                          element)
        self._add_element("fluxDensityLimit", metrics.flux_density_limit,
                          element)
        self._add_element("magLimit", metrics.mag_limit, element)
        if self._output_version >= 24:
            self._add_element("sampleSNR", metrics.sample_snr, element)

    def _add_quality_element(self, quality, parent):
        if quality is None:
            return

        element = self._get_caom_element("quality", parent)
        self._add_element("flag", quality.flag.value, element)

    def _add_observable_element(self, observable, parent):
        if observable is None:
            return

        element = self._get_caom_element("observable", parent)
        self._add_element("ucd", observable.ucd, element)
        if observable.calibration:
            if self._output_version < 25:
                raise AttributeError(
                    "Attempt to write CAOM2.5 element (observable.calibration) as "
                    "{} Observation".format(self._output_version))
            else:
                self._add_element("calibration", observable.calibration.value, element)

    def _add_transition_element(self, transition, parent):
        if transition is None:
            return

        element = self._get_caom_element("transition", parent)
        self._add_element("species", transition.species, element)
        self._add_element("transition", transition.transition, element)

    def _add_artifacts_element(self, artifacts, parent):
        if artifacts is None:
            return

        if self._output_format == 'xml':
            element = self._get_caom_element("artifacts", parent)
        else:
            element = parent["artifacts"] = []
        for _artifact in artifacts.values():
            if self._output_format == 'xml':
                artifact_element = self._get_caom_element("artifact", element)
            else:
                artifact_element = {}
                element.append(artifact_element)
            self._add_entity_attributes(_artifact, artifact_element)
            self._add_element("uri", _artifact.uri, artifact_element)
            if self._output_version >= 25:
                self._add_element("uriBucket", _artifact.uri_bucket, artifact_element)
            self._add_element("productType", _artifact.product_type.value,
                              artifact_element)
            self._add_element("releaseType", _artifact.release_type.value,
                              artifact_element)
            if self._output_version < 24 and _artifact.content_release:
                raise AttributeError(
                    "Attempt to write CAOM2.4 element "
                    "(artifact.content_realease) as CAOM2.3 Observation")
            self._add_datetime_element('contentRelease',
                                       _artifact.content_release,
                                       artifact_element)
            if self._output_version < 24 and _artifact.content_read_groups:
                raise AttributeError(
                    "Attempt to write CAOM2.4 element "
                    "(artifact.content_read_groups) as CAOM2.3 Observation")
            self._add_groups_element("contentReadGroups",
                                     _artifact.content_read_groups,
                                     artifact_element)
            self._add_element("contentType", _artifact.content_type,
                              artifact_element)
            self._add_element("contentLength", _artifact.content_length,
                              artifact_element)
            if _artifact.content_checksum:
                self._add_element("contentChecksum",
                                  _artifact.content_checksum,
                                  artifact_element)
            if _artifact.description_id is not None:
                if self._output_version < 25:
                    raise AttributeError(
                        "Attempt to write CAOM2.5 element "
                        "(artifact.description_id) as CAOM{} Observation".format(self._output_version))
                else:
                    self._add_element("descriptionID",
                                      _artifact.description_id, artifact_element)
            self._add_parts_element(_artifact.parts, artifact_element)

    def _add_parts_element(self, parts, parent):
        if not parts:
            return

        element = self._get_caom_element("parts", parent)
        for _part in parts.values():
            part_element = self._get_caom_element("part", element)
            self._add_entity_attributes(_part, part_element)
            self._add_element("name", _part.name, part_element)
            if _part.product_type is not None:
                self._add_element("productType", _part.product_type.value,
                                  part_element)
            self._add_chunks_element(_part.chunks, part_element)

    def _add_chunks_element(self, chunks, parent):
        if not chunks:
            return

        element = self._get_caom_element("chunks", parent)
        for _chunk in chunks:
            chunk_element = self._get_caom_element("chunk", element)
            self._add_entity_attributes(_chunk, chunk_element)
            if _chunk.product_type is not None:
                self._add_element("productType",
                                  _chunk.product_type.value,
                                  chunk_element)
            self._add_element("naxis", _chunk.naxis, chunk_element)
            self._add_element("observableAxis", _chunk.observable_axis,
                              chunk_element)
            self._add_element("positionAxis1", _chunk.position_axis_1,
                              chunk_element)
            self._add_element("positionAxis2", _chunk.position_axis_2,
                              chunk_element)
            self._add_element("energyAxis", _chunk.energy_axis, chunk_element)
            self._add_element("timeAxis", _chunk.time_axis, chunk_element)
            self._add_element("polarizationAxis", _chunk.polarization_axis,
                              chunk_element)
            if _chunk.custom_axis is not None:
                if self._output_version < 24:
                    raise AttributeError(
                        'Chunk.custom_axis only supported in CAOM2.4')
                else:
                    self._add_element("customAxis", _chunk.custom_axis,
                                      chunk_element)
            self._add_observable_axis_element(_chunk.observable, chunk_element)
            self._add_spatial_wcs_element(_chunk.position, chunk_element)
            self._add_spectral_wcs_element(_chunk.energy, chunk_element)
            self._add_temporal_wcs_element(_chunk.time, chunk_element)
            self._add_polarization_wcs_element(_chunk.polarization,
                                               chunk_element)
            if _chunk.custom is not None:
                if self._output_version < 24:
                    raise AttributeError(
                        'Chunk.custom only supported in CAOM2.4')
                else:
                    self._add_custom_wcs_element(_chunk.custom,
                                                 chunk_element)

    def _add_observable_axis_element(self, observable, parent):
        if observable is None:
            return

        element = self._get_caom_element("observable", parent)
        self._add_slice_element("dependent", observable.dependent, element)
        self._add_slice_element("independent", observable.independent, element)

    def _add_spatial_wcs_element(self, position, parent):
        """ Builds a representation of a SpatialWCS and adds it to the
            parent element. """
        if position is None:
            return

        element = self._get_caom_element("position", parent)
        self._add_coord_axis2d_element("axis", position.axis, element)
        self._add_element("coordsys", position.coordsys, element)
        self._add_element("equinox", position.equinox, element)
        self._add_element("resolution", position.resolution, element)

    def _add_spectral_wcs_element(self, energy, parent):
        """ Builds a representation of a SpectralWCS and adds it to the
            parent element."""
        if energy is None:
            return

        element = self._get_caom_element("energy", parent)
        self._add_coord_axis1d_element("axis", energy.axis, element)
        self._add_element("specsys", energy.specsys, element)
        self._add_element("ssysobs", energy.ssysobs, element)
        self._add_element("ssyssrc", energy.ssyssrc, element)
        self._add_element("restfrq", energy.restfrq, element)
        self._add_element("restwav", energy.restwav, element)
        self._add_element("velosys", energy.velosys, element)
        self._add_element("zsource", energy.zsource, element)
        self._add_element("velang", energy.velang, element)
        self._add_element("bandpassName", energy.bandpass_name, element)
        self._add_element("resolvingPower", energy.resolving_power, element)
        self._add_transition_element(energy.transition, element)

    def _add_temporal_wcs_element(self, time, parent):
        """ Builds a representation of a TemporalWCS and adds it to the
            parent element. """
        if time is None:
            return

        element = self._get_caom_element("time", parent)
        self._add_coord_axis1d_element("axis", time.axis, element)
        self._add_element("timesys", time.timesys, element)
        self._add_element("trefpos", time.trefpos, element)
        self._add_element("mjdref", time.mjdref, element)
        self._add_element("exposure", time.exposure, element)
        self._add_element("resolution", time.resolution, element)

    def _add_polarization_wcs_element(self, polarization, parent):
        """ Builds a representation of a PolarizationWCS and adds it to the
            parent element. """
        if polarization is None:
            return

        element = self._get_caom_element("polarization", parent)
        self._add_coord_axis1d_element("axis", polarization.axis, element)

    def _add_custom_wcs_element(self, custom, parent):
        """ Builds a representation of a CustomWCS and adds it to the
            parent element. """
        if custom is None:
            return

        element = self._get_caom_element("custom", parent)
        self._add_coord_axis1d_element("axis", custom.axis, element)

    # /*+ CAOM2 Types #-*/

    def _add_point_element(self, name, point, parent):
        """ Builds a representation of a Point and adds it to the
            parent element. """
        if point is None:
            return

        if self._output_format == 'xml':
            element = self._get_caom_element(name, parent)
        else:
            element = {}
            if isinstance(parent, list):
                parent.append(element)
            else:
                parent[name] = element

        self._add_element("cval1", point.cval1, element)
        self._add_element("cval2", point.cval2, element)

    # /*+ WCS Types #-*/

    def _add_axis_element(self, name, axis, parent):
        """ Builds a representation of a Axis and adds it to the
            parent element. """
        if axis is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_element("ctype", axis.ctype, element)
        if axis.cunit:
            self._add_element("cunit", axis.cunit, element)

    def _add_coord2d_element(self, name, coord, parent):
        """ Builds a representation of a Coord2D and adds it to the
            parent element. """
        if coord is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_ref_coord_element("coord1", coord.coord1, element)
        self._add_ref_coord_element("coord2", coord.coord2, element)

    def _add_value_coord2d_element(self, name, coord, parent):
        """ Builds a representation of a ValueCoord2D and adds it to the
            parent element. """
        if coord is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_element("coord1", coord.coord1, element)
        self._add_element("coord2", coord.coord2, element)

    def _add_coord_axis1d_element(self, name, axis, parent):
        """ Builds a representation of a CoordAxis1D and adds it to the
            parent element. """
        if axis is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_axis_element("axis", axis.axis, element)
        self._add_coord_error_element("error", axis.error, element)
        self._add_coord_range1d_element("range", axis.range, element)
        self._add_coord_bounds1d_element("bounds", axis.bounds, element)
        self._add_coord_function1d_element("function", axis.function, element)

    def _add_coord_axis2d_element(self, name, axis, parent):
        """ Builds a representation of a CoordAxis2D and adds it to the
            parent element. """
        if axis is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_axis_element("axis1", axis.axis1, element)
        self._add_axis_element("axis2", axis.axis2, element)
        self._add_coord_error_element("error1", axis.error1, element)
        self._add_coord_error_element("error2", axis.error2, element)
        self._add_coord_range2d_element("range", axis.range, element)
        self._add_coord_bounds2d_element("bounds", axis.bounds, element)
        self._add_coord_function2d_element("function", axis.function, element)

    def _add_coord_bounds1d_element(self, name, bounds, parent):
        """ Builds a representation of a CoordBounds1D and adds it to the
            parent element. """
        if bounds is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_coord_range_1d_list_element("samples", bounds.samples,
                                              element)

    def _add_coord_bounds2d_element(self, name, bounds, parent):
        """Builds a representation of a CoordBounds2D and adds it to the
            parent element. """
        if bounds is None:
            return

        element = self._get_caom_element(name, parent)
        if isinstance(bounds, wcs.CoordCircle2D):
            self._add_coord_circle2d_element("circle",
                                             wcs.CoordCircle2D(bounds.center,
                                                               bounds.radius),
                                             element)
        elif isinstance(bounds, wcs.CoordPolygon2D):
            self._add_coord_polygon2d_element("polygon", bounds, element)
        else:
            raise TypeError("BUG: unsupported CoordBounds2D type "
                            + bounds.__class__.__name__)

    def _add_coord_circle2d_element(self, name, circle, parent):
        """ Builds a representation of a CoordCircle2D and adds it to the
            parent element. """
        if circle is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_value_coord2d_element("center", circle.center, element)
        self._add_element("radius", circle.radius, element)

    def _add_coord_error_element(self, name, error, parent):
        """ Builds a representation of a CoordError and adds it to the
            parent element. """
        if error is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_element("syser", error.syser, element)
        self._add_element("rnder", error.rnder, element)

    def _add_coord_function1d_element(self, name, function, parent):
        """ Builds a representation of a CoordFunction1D and adds it to the
            parent element. """
        if function is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_element("naxis", function.naxis, element)
        self._add_element("delta", function.delta, element)
        self._add_ref_coord_element("refCoord", function.ref_coord, element)

    def _add_coord_function2d_element(self, name, function, parent):
        """ Builds a representation of a CoordFunction2D and adds it to the
            parent element. """
        if function is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_dimension2d_element("dimension", function.dimension, element)
        self._add_coord2d_element("refCoord", function.ref_coord, element)
        self._add_element("cd11", function.cd11, element)
        self._add_element("cd12", function.cd12, element)
        self._add_element("cd21", function.cd21, element)
        self._add_element("cd22", function.cd22, element)

    def _add_coord_polygon2d_element(self, name, polygon, parent):
        """ Builds a representation of a CoordPolygon2D and adds it to the
            parent element. """
        if polygon is None:
            return

        element = self._get_caom_element(name, parent)
        if len(polygon.vertices) > 0:
            vertices_element = self._get_caom_element("vertices", element)
            for vertex in polygon.vertices:
                self._add_value_coord2d_element("vertex", vertex,
                                                vertices_element)

    def _add_coord_range1d_element(self, name, _range, parent):
        """ Builds a representation of a CoordRange1D and adds it to the
            parent element. """
        if _range is None:
            return
        element = self._get_caom_element(name, parent)
        self._add_ref_coord_element("start", _range.start, element)
        self._add_ref_coord_element("end", _range.end, element)

    def _add_coord_range2d_element(self, name, _range, parent):
        """ Builds a representation of a CoordRange2D and adds it to the
            parent element. """
        if _range is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_coord2d_element("start", _range.start, element)
        self._add_coord2d_element("end", _range.end, element)

    def _add_dimension2d_element(self, name, dimension, parent):
        """ Builds a representation of a Dimension2D and adds it to the
            parent element. """
        if dimension is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_element("naxis1", dimension.naxis1, element)
        self._add_element("naxis2", dimension.naxis2, element)

    def _add_ref_coord_element(self, name, ref_coord, parent):
        """ Builds a representation of a RefCoord and adds it to the
            parent element. """
        if ref_coord is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_element("pix", ref_coord.pix, element)
        self._add_element("val", ref_coord.val, element)

    def _add_slice_element(self, name, _slice, parent):
        """ Builds a representation of a Slice and adds it to the
            parent element. """
        if _slice is None:
            return

        element = self._get_caom_element(name, parent)
        self._add_axis_element("axis", _slice.axis, element)
        self._add_element("bin", _slice.bin, element)

    def _add_attribute(self, name, value, element):
        if self._output_format == 'xml':
            element.set(self._caom2_namespace + name, value)
        else:
            element['@' + name] = value

    def _add_element(self, name, value, parent):
        if value is None:
            return
        if self._output_format == 'xml':
            element = self._get_caom_element(name, parent)
            if isinstance(value, str):
                element.text = value
            elif isinstance(value, Enum):
                element.text = value.value
            else:
                element.text = str(value)
        else:
            if not isinstance(value, (str, int, float, bool)):
                raise TypeError("CAOM JSON value must be str, int, float or bool not "
                                + value.__class__.__name__)
            if isinstance(parent, list):
                parent.append(value)
            else:
                parent[name] = value

    def _add_boolean_element(self, name, value, parent):
        if value is None:
            return
        if self._output_format == 'xml':
            element = self._get_caom_element(name, parent)
            element.text = str(value).lower()
        else:
            parent[name] = bool(value)

    def _add_datetime_element(self, name, value, parent):
        if value is None:
            return
        if self._output_format == 'xml':
            element = self._get_caom_element(name, parent)
            element.text = caom_util.date2ivoa(value)
        else:
            parent[name] = caom_util.date2ivoa(value)

    def _add_keywords_element(self, collection, parent):
        if collection is None or \
                (len(collection) == 0 and not self._write_empty_collections):
            return
        if self._output_format == 'xml':
            element = self._get_caom_element("keywords", parent)
            for keyword in collection:
                self._get_caom_element("keyword", element).text = keyword
        else:
            parent['keywords'] = sorted(list(collection))

    def _add_coord_range_1d_list_element(self, name, values, parent):
        if values is None:
            return
        if self._output_format == 'xml':
            element = self._get_caom_element(name, parent)
            for v in values:
                self._add_coord_range1d_element("range", v, element)
        else:
            raise NotImplementedError(
                "CAOM JSON output not implemented")

    def _add_inputs_element(self, name, collection, parent):
        if collection is None or \
                (len(collection) == 0 and not self._write_empty_collections):
            return
        element = self._get_caom_element(name, parent)
        if self._output_version < 25:
            input_tag = "planeURI"
        else:
            input_tag = "input"
        if self._output_format == 'xml':
            for plane_uri in collection:
                self._add_element(input_tag, plane_uri, element)
        else:
            parent[name] = [inp for inp in collection]

    def _get_caom_element(self, tag, parent):
        if self._output_format == 'xml':
            return etree.SubElement(parent, self._caom2_namespace + tag)
        else:
            parent[tag] = {}
            return parent[tag]


class ObservationParsingException(Exception):
    pass
