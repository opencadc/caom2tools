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

""" Defines ObservationWriter class """

from lxml import etree
import pkg_resources
from caom2_observation_reader import ObservationReader
from .. caom2_observation import Observation
from .. caom2_composite_observation import CompositeObservation
from .. caom2_simple_observation import SimpleObservation
from .. wcs.caom2_coord_circle2d import CoordCircle2D
from .. wcs.caom2_coord_polygon2d import CoordPolygon2D
from .. caom2_enums import CalibrationLevel
from .. caom2_enums import DataProductType
from .. caom2_enums import ObservationIntentType
from .. caom2_enums import ProductType
from .. caom2_enums import TargetType
from .. util.caom2_util import date2ivoa


class ObservationWriter(object):
    """ ObservationWriter """

    CAOM2_NAMESPACE = "vos://cadc.nrc.ca!vospace/CADC/xml/CAOM/v2.0"
    CAOM2 = "{%s}" % CAOM2_NAMESPACE
    XSI_NAMESPACE = "http://www.w3.org/2001/XMLSchema-instance"
    XSI = "{%s}" % XSI_NAMESPACE
    NSMAP = {"caom2": CAOM2_NAMESPACE, "xsi": XSI_NAMESPACE}

    def __init__(self, validate=False, write_empty_collections=False):
        self._validate = validate
        self._write_empty_collections = write_empty_collections

        schema_path = pkg_resources.resource_filename(
            ObservationReader.CAOM2_PKG, ObservationReader.SCHEMA_FILE)
        xmlschema_doc = etree.parse(schema_path)
        self._xmlschema = etree.XMLSchema(xmlschema_doc)

    def write(self, observation, out):
        assert isinstance(observation, Observation), (
            "observation is not an Observation")

        obs = etree.Element(self.CAOM2 + "Observation", nsmap=self.NSMAP)
        if isinstance(observation, SimpleObservation):
            obs.set(self.XSI + "type", "caom2:SimpleObservation")
        else:
            obs.set(self.XSI + "type", "caom2:CompositeObservation")

        self._addAttribute("id", str(observation._id), obs)
        if (observation._last_modified != None):
            self._addAttribute("lastModified",
                               date2ivoa(observation._last_modified), obs)

        self._addElement("collection", observation.collection, obs)
        self._addElement("observationID", observation.observation_id, obs)
        self._addDatetimeElement("metaRelease", observation.meta_release, obs)
        self._addElement("sequenceNumber", observation.sequence_number, obs)
        self._addAlgorithmElement(observation.algorithm, obs)
        self._addElement("type", observation.obs_type, obs)
        if (observation.intent != None):
            self._addElement("intent",
                ObservationIntentType.get(str(observation.intent)).value, obs)

        self._addProposalElement(observation.proposal, obs)
        self._addTargetElement(observation.target, obs)
        self._addTelescopeElement(observation.telescope, obs)
        self._addInstrumentElement(observation.instrument, obs)
        self._addEnvironmentElement(observation.environment, obs)
        self._addPlanesElement(observation.planes, obs)

        if isinstance(observation, CompositeObservation):
            self._addMembersElement(observation.members, obs)

        if (self._validate):
            self._xmlschema.validate(obs)

        out.write(etree.tostring(obs, xml_declaration=True, encoding='UTF-8',
                                 pretty_print=True))

    def _addAlgorithmElement(self, algorithm, parent):
        if (algorithm is None):
            return

        element = self._getCaom2Element("algorithm", parent)
        self._addElement("name", algorithm.name, element)

    def _addProposalElement(self, proposal, parent):
        if (proposal == None):
            return

        element = self._getCaom2Element("proposal", parent)
        self._addElement("id", proposal.proposal_id, element)
        self._addElement("pi", proposal.pi_name, element)
        self._addElement("project", proposal.project, element)
        self._addElement("title", proposal.title, element)
        self._addListElement("keywords", proposal.keywords, element)

    def _addTargetElement(self, target, parent):
        if (target == None):
            return

        element = self._getCaom2Element("target", parent)
        self._addElement("name", target.name, element)
        if (target.target_type != None):
            self._addElement("type",
                TargetType.get(str(target.target_type)).value, element)
        if (target.standard != None):
            self._addElement("standard", str(target.standard).lower(), element)
        self._addElement("redshift", target.redshift, element)
        self._addListElement("keywords", target.keywords, element)

    def _addTelescopeElement(self, telescope, parent):
        if (telescope == None):
            return

        element = self._getCaom2Element("telescope", parent)
        self._addElement("name", telescope.name, element)
        self._addElement("geoLocationX", telescope.geo_location_x, element)
        self._addElement("geoLocationY", telescope.geo_location_y, element)
        self._addElement("geoLocationZ", telescope.geo_location_z, element)
        self._addListElement("keywords", telescope.keywords, element)

    def _addInstrumentElement(self, instrument, parent):
        if (instrument == None):
            return

        element = self._getCaom2Element("instrument", parent)
        self._addElement("name", instrument.name, element)
        self._addListElement("keywords", instrument.keywords, element)

    def _addEnvironmentElement(self, environment, parent):
        if (environment == None):
            return

        element = self._getCaom2Element("environment", parent)
        self._addElement("seeing", environment.seeing, element)
        self._addElement("humidity", environment.humidity, element)
        self._addElement("elevation", environment.elevation, element)
        self._addElement("tau", environment.tau, element)
        self._addElement("wavelengthTau", environment.wavelength_tau, element)
        self._addElement("ambientTemp", environment.ambient_temp, element)
        if (environment.photometric != None):
            self._addElement("photometric",
                             str(environment.photometric).lower(), element)

    def _addMembersElement(self, members, parent):
        if (members == None or
            (len(members) == 0 and not self._write_empty_collections)):
            return

        element = self._getCaom2Element("members", parent)
        for member in members:
            memberElement = self._getCaom2Element("observationURI", element)
            memberElement.text = member.uri

    def _addPlanesElement(self, planes, parent):
        if (planes == None or
            (len(planes) == 0 and not self._write_empty_collections)):
            return

        element = self._getCaom2Element("planes", parent)
        for plane in planes.itervalues():
            planeElement = self._getCaom2Element("plane", element)
            self._addAttribute("id", str(plane._id), planeElement)
            if (plane._last_modified != None):
                self._addAttribute("lastModified",
                              date2ivoa(plane._last_modified), planeElement)
            self._addElement("productID", plane.product_id, planeElement)
            self._addDatetimeElement("metaRelease", plane.meta_release,
                                    planeElement)
            self._addDatetimeElement("dataRelease", plane.data_release,
                                    planeElement)
            if (plane.data_product_type != None):
                self._addElement("dataProductType",
                    DataProductType.get(str(plane.data_product_type)).value,
                    planeElement)
            if (plane.calibration_level != None):
                self._addElement("calibrationLevel",
                    CalibrationLevel(str(plane.calibration_level)).value,
                    planeElement)
            self._addProvenanceElement(plane.provenance, planeElement)
            self._addMetricsElement(plane.metrics, planeElement)
            self._addArtifactsElement(plane.artifacts, planeElement)

    def _addProvenanceElement(self, provenance, parent):
        if (provenance == None):
            return

        element = self._getCaom2Element("provenance", parent)
        self._addElement("name", provenance.name, element)
        self._addElement("version", provenance.version, element)
        self._addElement("project", provenance.project, element)
        self._addElement("producer", provenance.producer, element)
        self._addElement("runID", provenance.run_id, element)
        self._addElement("reference", provenance.reference, element)
        self._addDatetimeElement("lastExecuted", provenance.last_executed,
                                element)
        self._addListElement("keywords", provenance.keywords, element)
        self._addInputsElement("inputs", provenance.inputs, element)

    def _addMetricsElement(self, metrics, parent):
        if (metrics == None):
            return

        element = self._getCaom2Element("metrics", parent)
        self._addElement("sourceNumberDensity", metrics.source_number_density,
                        element)
        self._addElement("background", metrics.background, element)
        self._addElement("backgroundStddev", metrics.background_std_dev,
                        element)
        self._addElement("fluxDensityLimit", metrics.flux_density_limit,
                        element)
        self._addElement("magLimit", metrics.mag_limit, element)

    def _addTransitionElement(self, transition, parent):
        if (transition == None):
            return

        element = self._getCaom2Element("transition", parent)
        self._addElement("species", transition.species, element)
        self._addElement("transition", transition.transition, element)

    def _addArtifactsElement(self, artifacts, parent):
        if (artifacts == None):
            return

        element = self._getCaom2Element("artifacts", parent)
        for artifact in artifacts.itervalues():
            artifactElement = self._getCaom2Element("artifact", element)
            self._addAttribute("id", str(artifact._id), artifactElement)
            if (artifact._last_modified != None):
                self._addAttribute("lastModified",
                              date2ivoa(artifact._last_modified),
                              artifactElement)
            self._addElement("uri", artifact.uri, artifactElement)
            self._addElement("contentType", artifact.content_type,
                            artifactElement)
            self._addElement("contentLength", artifact.content_length,
                            artifactElement)
            if (artifact.product_type != None):
                self._addElement("productType",
                    ProductType.get(str(artifact.product_type)).value,
                    artifactElement)
            self._addElement("alternative", str(artifact.alternative).lower(),
                            artifactElement)
            self._addPartsElement(artifact.parts, artifactElement)

    def _addPartsElement(self, parts, parent):
        if (parts == None):
            return

        element = self._getCaom2Element("parts", parent)
        for part in parts.itervalues():
            partElement = self._getCaom2Element("part", element)
            self._addAttribute("id", str(part._id), partElement)
            if (part._last_modified != None):
                self._addAttribute("lastModified",
                                   date2ivoa(part._last_modified),
                                   partElement)
            self._addElement("name", part.name, partElement)
            if (part.product_type != None):
                self._addElement("productType",
                    ProductType.get(str(part.product_type)).value, partElement)
            self._addChunksElement(part.chunks, partElement)

    def _addChunksElement(self, chunks, parent):
        if (chunks == None):
            return

        element = self._getCaom2Element("chunks", parent)
        for chunk in chunks:
            chunkElement = self._getCaom2Element("chunk", element)
            self._addAttribute("id", str(chunk._id), chunkElement)
            if (chunk._last_modified != None):
                self._addAttribute("lastModified",
                                   date2ivoa(chunk._last_modified),
                                   chunkElement)
            if (chunk.product_type != None):
                self._addElement("productType",
                    ProductType.get(str(chunk.product_type)).value,
                    chunkElement)
            self._addElement("naxis", chunk.naxis, chunkElement)
            self._addElement("observableAxis", chunk.observable_axis,
                            chunkElement)
            self._addElement("positionAxis1", chunk.position_axis_1,
                            chunkElement)
            self._addElement("positionAxis2", chunk.position_axis_2,
                            chunkElement)
            self._addElement("energyAxis", chunk.energy_axis, chunkElement)
            self._addElement("timeAxis", chunk.time_axis, chunkElement)
            self._addElement("polarizationAxis", chunk.polarization_axis,
                            chunkElement)

            self._addObservableAxisElement(chunk.observable, chunkElement)
            self._addSpatialWCSElement(chunk.position, chunkElement)
            self._addSpectralWCSElement(chunk.energy, chunkElement)
            self._addTemporalWCSElement(chunk.time, chunkElement)
            self._addPolarizationWCSElement(chunk.polarization, chunkElement)

    def _addObservableAxisElement(self, observable, parent):
        if (observable == None):
            return

        element = self._getCaom2Element("observable", parent)
        self._addSliceElement("dependent", observable.dependent, element)
        self._addSliceElement("independent", observable.independent, element)

    def _addSpatialWCSElement(self, position, parent):
        """ Builds a representation of a SpatialWCS and adds it to the
            parent element. """
        if (position == None):
            return

        element = self._getCaom2Element("position", parent)
        self._addCoordAxis2DElement("axis", position.axis, element)
        self._addElement("coordsys", position.coordsys, element)
        self._addElement("equinox", position.equinox, element)
        self._addElement("resolution", position.resolution, element)

    def _addSpectralWCSElement(self, energy, parent):
        """ Builds a representation of a SpectralWCS and adds it to the
            parent element."""
        if (energy == None):
            return

        element = self._getCaom2Element("energy", parent)
        self._addCoordAxis1DElement("axis", energy.axis, element)
        self._addElement("specsys", energy.specsys, element)
        self._addElement("ssysobs", energy.ssysobs, element)
        self._addElement("ssyssrc", energy.ssyssrc, element)
        self._addElement("restfrq", energy.restfrq, element)
        self._addElement("restwav", energy.restwav, element)
        self._addElement("velosys", energy.velosys, element)
        self._addElement("zsource", energy.zsource, element)
        self._addElement("velang", energy.velang, element)
        self._addElement("bandpassName", energy.bandpass_name, element)
        self._addElement("resolvingPower", energy.resolving_power, element)
        self._addTransitionElement(energy.transition, element)

    def _addTemporalWCSElement(self, time, parent):
        """ Builds a representation of a TemporalWCS and adds it to the
            parent element. """
        if (time == None):
            return

        element = self._getCaom2Element("time", parent)
        self._addCoordAxis1DElement("axis", time.axis, element)
        self._addElement("timesys", time.timesys, element)
        self._addElement("trefpos", time.trefpos, element)
        self._addElement("mjdref", time.mjdref, element)
        self._addElement("exposure", time.exposure, element)
        self._addElement("resolution", time.resolution, element)

    def _addPolarizationWCSElement(self, polarization, parent):
        """ Builds a representation of a PolarizationWCS and adds it to the
            parent element. """
        if (polarization == None):
            return

        element = self._getCaom2Element("polarization", parent)
        self._addCoordAxis1DElement("axis", polarization.axis, element)

#/*+ WCS Types #-*/

    def _addAxisElement(self, name, axis, parent):
        """ Builds a representation of a Axis and adds it to the
            parent element. """
        if (axis == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addElement("ctype", axis.ctype, element)
        if (axis.cunit):
            self._addElement("cunit", axis.cunit, element)

    def _addCoord2DElement(self, name, coord, parent):
        """ Builds a representation of a Coord2D and adds it to the
            parent element. """
        if(coord == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addRefCoordElement("coord1", coord.coord1, element)
        self._addRefCoordElement("coord2", coord.coord2, element)

    def _addCoordAxis1DElement(self, name, axis, parent):
        """ Builds a representation of a CoordAxis1D and adds it to the
            parent element. """
        if (axis == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addAxisElement("axis", axis.axis, element)
        self._addCoordErrorElement("error", axis.error, element)
        self._addCoordRange1DElement("range", axis.range, element)
        self._addCoordBounds1DElement("bounds", axis.bounds, element)
        self._addCoordFunction1DElement("function", axis.function, element)

    def _addCoordAxis2DElement(self, name, axis, parent):
        """ Builds a representation of a CoordAxis2D and adds it to the
            parent element. """
        if (axis == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addAxisElement("axis1", axis.axis1, element)
        self._addAxisElement("axis2", axis.axis2, element)
        self._addCoordErrorElement("error1", axis.error1, element)
        self._addCoordErrorElement("error2", axis.error2, element)
        self._addCoordRange2DElement("range", axis.range, element)
        self._addCoordBounds2DElement("bounds", axis.bounds, element)
        self._addCoordFunction2DElement("function", axis.function, element)

    def _addCoordBounds1DElement(self, name, bounds, parent):
        """ Builds a representation of a CoordBounds1D and adds it to the
            parent element. """
        if (bounds == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addCoordRange1DListElement("samples", bounds.samples, element)

    def _addCoordBounds2DElement(self, name, bounds, parent):
        """Builds a representation of a CoordBounds2D and adds it to the
            parent element. """
        if (bounds == None):
            return

        element = self._getCaom2Element(name, parent)
        if isinstance(bounds, CoordCircle2D):
            self._addCoordCircle2DElement("circle",
                                         CoordCircle2D(bounds.center,
                                                       bounds.radius), element)
        elif isinstance(bounds, CoordPolygon2D):
            self._addCoordPolygon2DElement("polygon", bounds, element)
        else:
            raise TypeError("BUG: unsupported CoordBounds2D type "
                            + bounds.__class__.__name__)

    def _addCoordCircle2DElement(self, name, circle, parent):
        """ Builds a representation of a CoordCircle2D and adds it to the
            parent element. """
        if (circle == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addCoord2DElement("center", circle.center, element)
        self._addElement("radius", circle.radius, element)

    def _addCoordErrorElement(self, name, error, parent):
        """ Builds a representation of a CoordError and adds it to the
            parent element. """
        if (error == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addElement("syser", error.syser, element)
        self._addElement("rnder", error.rnder, element)

    def _addCoordFunction1DElement(self, name, function, parent):
        """ Builds a representation of a CoordFunction1D and adds it to the
            parent element. """
        if (function == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addElement("naxis", function.naxis, element)
        self._addElement("delta", function.delta, element)
        self._addRefCoordElement("refCoord", function.ref_coord, element)

    def _addCoordFunction2DElement(self, name, function, parent):
        """ Builds a representation of a CoordFunction2D and adds it to the
            parent element. """
        if (function == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addDimension2DElement("dimension", function.dimension, element)
        self._addCoord2DElement("refCoord", function.ref_coord, element)
        self._addElement("cd11", function.cd11, element)
        self._addElement("cd12", function.cd12, element)
        self._addElement("cd21", function.cd21, element)
        self._addElement("cd22", function.cd22, element)

    def _addCoordPolygon2DElement(self, name, polygon, parent):
        """ Builds a representation of a CoordPolygon2D and adds it to the
            parent element. """
        if (polygon == None):
            return

        element = self._getCaom2Element(name, parent)
        if (len(polygon.vertices) > 0):
            verticesElement = self._getCaom2Element("vertices", element)
            for vertex in polygon.vertices:
                self._addCoord2DElement("vertex", vertex, verticesElement)

    def _addCoordRange1DElement(self, name, _range, parent):
        """ Builds a representation of a CoordRange1D and adds it to the
            parent element. """
        if (_range == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addRefCoordElement("start", _range.start, element)
        self._addRefCoordElement("end", _range.end, element)

    def _addCoordRange2DElement(self, name, _range, parent):
        """ Builds a representation of a CoordRange2D and adds it to the
            parent element. """
        if (_range == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addCoord2DElement("start", _range.start, element)
        self._addCoord2DElement("end", _range.end, element)

    def _addDimension2DElement(self, name, dimension, parent):
        """ Builds a representation of a Dimension2D and adds it to the
            parent element. """
        if (dimension == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addElement("naxis1", dimension.naxis1, element)
        self._addElement("naxis2", dimension.naxis2, element)

    def _addRefCoordElement(self, name, refCoord, parent):
        """ Builds a representation of a RefCoord and adds it to the
            parent element. """
        if (refCoord == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addElement("pix", refCoord.pix, element)
        self._addElement("val", refCoord.val, element)

    def _addSliceElement(self, name, _slice, parent):
        """ Builds a representation of a Slice and adds it to the
            parent element. """
        if (_slice == None):
            return

        element = self._getCaom2Element(name, parent)
        self._addAxisElement("axis", _slice.axis, element)
        self._addElement("bin", _slice.bin, element)

    def _addAttribute(self, name, value, element):
        element.set(self.CAOM2 + name, value)

    def _addElement(self, name, text, parent):
        if (text == None):
            return
        element = self._getCaom2Element(name, parent)
        if (isinstance(text, str)):
            element.text = text
        else:
            element.text = str(text)

    def _addDatetimeElement(self, name, value, parent):
        if (value == None):
            return
        element = self._getCaom2Element(name, parent)
        element.text = date2ivoa(value)

    def _addListElement(self, name, collection, parent):
        if (collection == None or
            (len(collection) == 0 and not self._write_empty_collections)):
            return
        element = self._getCaom2Element(name, parent)
        element.text = ' '.join(collection)

    def _addCoordRange1DListElement(self, name, values, parent):
        if (values == None):
            return
        element = self._getCaom2Element(name, parent)
        for v in values:
            self._addCoordRange1DElement("range", v, element)

    def _addInputsElement(self, name, collection, parent):
        if (collection == None or
            (len(collection) == 0 and not self._write_empty_collections)):
            return
        element = self._getCaom2Element(name, parent)
        for plane_uri in collection:
            self._addElement("planeURI", plane_uri.uri, element)

    def _getCaom2Element(self, tag, parent):
        return etree.SubElement(parent, self.CAOM2 + tag)
