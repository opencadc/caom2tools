This is the user documentation for the CADC Python Data Engineering Tools.

CAOM XML records can be created and modified using customized software built off the caom2 (TODO link) and caom2utils (TODO link) modules.
These modules enable the user to create more complex relations between observational datasets than can be achieved using the fits2caom2 (TODO link) process. The modules assist the user in creating CAOM objects which can then be written to an XML file and uploaded to the CADC caom2repo (TODO link) service.

# INSTALL

caom2 is a Python library, that can be installed using Python's pip command:

   ```
   pip install caom2
   ```

The source code is here (TODO link).

caom2utils is a Python library and the Python implementation of fits2caom2. It can be installed using Python's pip command:

   ```
   pip install caom2utils
   ```
   
The source code is here (TODO link).


# USAGE

caom2 (TODO link) and caom2utils (TODO link) are python modules that enable writing customized software. One must build the various objects that make up the Observation records, put them together, and then provide that object to the CAOM2 repository service (TODO link).

## Example

Here (TODO link to caom2 README.md) is a bare-bones example of how one might use caom2 (TODO link) to build a CAOM XML record.

Here (TODO link to opencadc-metadata-curation) are operational examples of workflows that incorporate caom2 (TODO link) and caom2utils (TODO link) to manage the entire life-cycle of a CAOM observation.
