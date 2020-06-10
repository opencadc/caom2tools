# Working With CAOM2

For observations to appear in [CADC search services](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/), an observation must first be described by a CAOM record. That description will then need to be loaded into the CADC CAOM repository, using a CADC web service. This web service will create a corresponding database record.

There are two ways to describe and load an Observation:
1. [Create a Python script](./script_description.md)
1. [Use the Command Line Interface](./cli_description.md)

Once an Observation has been described and loaded, it is searchable from CADC's UI.
