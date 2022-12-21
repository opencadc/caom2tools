caom2utils
==========

.. image:: https://img.shields.io/pypi/v/caom2utils.svg   
    :target: https://pypi.python.org/pypi/caom2utils

Utilities to facililate working with the CAOM2 model.

Observation Validation
----------------------

Validates a CAOM2 element (Observation, Plane, Artifact, Part or Chunk) with respect to the attributes of the element and possibly all its sub-elements. Example of validations: attribute values, spherical geometry of planes, WCS of chuncks, etc.

.. code:: python

    import sys
    from caom2 import SimpleObservation
    import caom2utils

    obs = SimpleObservation('collection', 'observationID')

    # change and update obs

    try:
        caom2utils.validate(obs)
    except Exception:
        print('My exception is not valid')


Observation Generation
======================

There are two command-line applications that will generate a CAOM2 Observation in this package. fits2caom2 assumes legacy operations, and should be a not-quite replacement for the Java application. caom2gen provides a bit different invocation mechanism, relying on blueprints for configuration information, and specifies the cardinality of the Observation/Plane/Artifact combinations at the command line.


fits2caom2
----------

.. code-block:: console

    usage: fits2caom2 [-h] [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                      [--host HOST] [--resource-id RESOURCE_ID] [-d | -q | -v]
                      [-V] [--dumpconfig] [--no_validate]
                      [-o OUT_OBS_XML]
                      (-i IN_OBS_XML | --observation collection observationID)
                      [--local LOCAL [LOCAL ...]] [--keep] [--test]
                      [--productID PRODUCTID] [--config CONFIG]
                      [--default DEFAULT] [--override OVERRIDE]
                      fileURI [fileURI ...]

    Augments an observation with information in one or more fits files.

    positional arguments:
      fileURI                                 URI of a fits file

    optional arguments:
      --cert CERT                             location of your X509 certificate to use for
                                              authentication (unencrypted, in PEM format)
      --config CONFIG                         optional CAOM2 utype to keyword config file to merge
                                              with the internal configuration
      -d, --debug                             debug messages
      --default DEFAULT                       file with default values for keywords
      --dumpconfig                            output the utype to keyword mapping to the console
      -h, --help                              show this help message and exit
      --host HOST                             base hostname for services - used mainly for testing
                                              (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
      -i, --in IN_OBS_XML                     input of observation to be augmented in XML
      --keep                                  keep the locally stored files after ingestion
      --local LOCAL [ LOCAL ...]
                                              list of files in local filesystem (same order as uri)
      -n                                      use .netrc in $HOME for authentication
      --netrc-file NETRC_FILE                 
                                              netrc file to use for authentication
      --no_validate                           by default, the application will validate the WCS
                                              information for an observation. Specifying this flag
                                              skips that step.
      --observation collection observationID  
                                              observation in a collection
      -o, --out OUT_OBS_XML                   output of augmented observation in XML
      --override OVERRIDE                     file with override values for keywords
      --productID PRODUCTID                   product ID of the plane in the observation
      -q, --quiet                             run quietly
      --resource-id RESOURCE_ID               resource identifier (default
                                              ivo://cadc.nrc.ca/fits2caom2)
      --test                                  test mode, do not persist to database
      -u, --user USER                         name of user to authenticate. Note: application
                                              prompts for the corresponding password!
      -v, --verbose                           verbose messages
      -V, --version                           show program's version number and exit


caom2gen
--------

.. code-block:: console

    usage: caom2gen [-h] [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER] [--host HOST] 
                [--resource-id RESOURCE_ID] [-d | -q | -v] [-V]
                [--dumpconfig] [--not_connected] [--no_validate] [-o OUT_OBS_XML] [--caom_namespace CAOM_NAMESPACE]
                (-i IN_OBS_XML | --observation collection observationID) [--local LOCAL [LOCAL ...]]
                [--external_url EXTERNAL_URL [EXTERNAL_URL ...]] [--module MODULE] [--plugin PLUGIN] 
                [--lineage LINEAGE [LINEAGE ...]]
                [--use_blueprint_parser USE_BLUEPRINT_PARSER [USE_BLUEPRINT_PARSER ...]] --blueprint BLUEPRINT [BLUEPRINT ...]

    Augments an observation with information in one or more fits files.

    optional arguments:
       --blueprint BLUEPRINT [BLUEPRINT ...]  list of files with blueprints for CAOM2 
                                              construction, in serialized format. If 
                                              the list is of length 1, the same 
                                              blueprint will be applied to all lineage 
                                              entries. Otherwise, there must be a 
                                              blueprint file per lineage entry.
       --caom_namespace CAOM_NAMESPACE        if this parameter is specified, over-ride 
                                              the default CAOM2 version when writing 
                                              XML. The default is the latest version of 
                                              CAOM2.3.
       --cert CERT                            location of your X509 certificate to use 
                                              for authentication (unencrypted, in PEM 
                                              format)
       -d, --debug                            debug messages
       --dumpconfig                           output the utype to keyword mapping to 
                                              the console
       --external_url EXTERNAL_URL [EXTERNAL_URL ...]
                                              service endpoint(s) that return(s) a 
                                              string that can be made into FITS 
                                              headers. Cardinality should be consistent 
                                              with lineage.
       -h, --help                             show this help message and exit
       --host HOST                            base hostname for services - used mainly 
                                              for testing (default: 
                                              www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
       -i, --in IN_OBS_XML                    input of observation to be augmented in XML
       --lineage LINEAGE [LINEAGE ...]        productID/artifactURI. List of plane/artifact 
                                              identifiers that will be created for the 
                                              identified observation.
       --local LOCAL [LOCAL ...]              list of files in local filesystem (same 
                                              order as uri)
       --module MODULE                        if the blueprint contains function calls, 
                                              call importlib.import_module for the named 
                                              module. Provide a fully qualified name. 
                                              Parameter choices are the artifact URI (uri) 
                                              or a list of astropy Header instances 
                                              (header). This will allow the update of a 
                                              single blueprint entry with a single call.
       -n                                     use .netrc in $HOME for authentication
       --netrc-file NETRC_FILE                netrc file to use for authentication
       --no_validate                          by default, the application will validate 
                                              the WCS information for an observation. 
                                              Specifying this flag skips that step.
       --not_connected                        if set, there is no internet connection, 
       so skip service invocations.
       --observation collection observationID observation in a collection
       -o, --out OUT_OBS_XML                  output of augmented observation in XML
       --plugin PLUGIN                        if this parameter is specified, call 
                                              importlib.import_module for the named 
                                              module. Then execute the method "update", 
                                              with the signature (Observation, **kwargs). 
                                              This will allow for the update of multiple 
                                              observation data members with one call.
       -q, --quiet                            run quietly
       --resource-id RESOURCE_ID              resource identifier (default 
                                              ivo://cadc.nrc.ca/fits2caom2)
       --use_blueprint_parser USE_BLUEPRINT_PARSER [USE_BLUEPRINT_PARSER ...]
                                              productID/artifactURI. List of lineage 
                                              entries that will be processed with a 
                                              BlueprintParser. Good for files with no
                                              metadata in the content.
       -u, --user USER                        name of user to authenticate. Note: 
                                              application prompts for the corresponding 
                                              password!
       -v, --verbose                          verbose messages
       -V, --version                          show program's version number and exit
