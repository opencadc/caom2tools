# Working With CAOM2

For observations to appear in [CADC search services](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/), an observation must first be described by a CAOM record. That description will then need to be loaded into the CADC CAOM repository, using a CADC web service. This web service will create a corresponding database record. 

Once an Observation has been described and loaded, it is searchable from CADC's UI.

* If you are interested in using CADC Python Data Engineering tools, you should start [here](./user/cli_description.md).
  
* If you are interested in scripting with the CADC Python Data Engineering tools, you should start [here](./user/script_description.md).
  
## Preconditions

1. These descriptions assume:
    1. a working knowledge of python. [Prefer python3, please](https://pythonclock.org/),
    1. a linux-type environment,
    1. a working directory location, where all files discussed are placed, and
    1. that you have a [CADC account](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/auth/request.html), which is configured by CADC to have read and write access to a CAOM `COLLECTION`.

1. This description uses the parameters `TEST_FILE.FITS`, `TEST_OBS.XML` and `COLLECTION`. Replace these values appropriately when executing the commands.

1. Copy the file `TEST_FILE.FITS` in the working directory. The metadata in this file will be described in the CAOM Observation created during this example.

1. The example will cause an instance to be created in the [CAOM2 sandbox](http://sc2.canfar.net/search/).  If you click the CAOM2 sandbox URL prior to the creation of the first CAOM instance for a `COLLECTION`, that `COLLECTION` will not show in the `Additional Constraints -> Collection` . Even after successful creation of a CAOM instance, it can take up to one day for the `COLLECTION` to be selectable from the UI. The CAOM2 'sandbox' is a site that mimics the production CADC CAOM2 storage service and search UI. This sandbox site allows developers and scientists to debug collection-specific code for creating and updating CAOM2 Observations. It also allows developers and scientiests to immediately view CAOM2 records as they will appear to users in the search interface.

    To use production CADC services, remove `resource-id` parameters in `caom2-repo` commands.

1. Install the following python dependencies:

    ```
    pip install caom2repo
    pip install caom2utils
    ```

1. Get credentials organized. The examples assume the use of a [./.netrc file](https://www.systutorials.com/docs/linux/man/5-netrc/). The examples expect this file to be named `./.netrc`, located in the working directory, and with permissions set to `-rw-------`. The `./.netrc` file content should include the following, with cadcusername and cadcpassword replaced with your CADC username and password values:

    ````
    machine www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca login canfarusername password canfarpassword
    machine www.canfar.net login canfarusername password canfarpassword
    machine ws-cadc.canfar.net login canfarusername password canfarpassword
    machine sc2.canfar.net login canfarusername password canfarpassword
    ````
    
    To set the `./.netrc` file permissions:
    
    ```
    chmod 600 ./.netrc
    ```

1. The caom2-repo client also supports username/password and X509 certificates. If you want to use X509 
certificates use the --cert parameter instead of the -n parameter in all the commands. The command line client `cadc-get-cert` is installed with the prerequisites for the `caom2repo` package, and `cadc-get-cert --help` from a terminal prompt will describe how to obtain a CADC certificate.


1. Test the install. Commands are case-sensitive.

    ```
    caom2-repo read --netrc ./.netrc --resource-id ivo://cadc.nrc.ca/sc2repo COLLECTION abc
    ```

    If the install was successful, this will report an error:

    ```
    Client Error: Not Found for url: http://sc2.canfar.net/sc2repo/auth-observations/COLLECTION/abc.
    ```

## Troubleshooting

1. If `pip install caom2utils` fails with the following error:

    ```
    AttributeError: module ‘enum’ has no attribute ‘IntFlag’
    ```

    Ensure the version of vos is >= 3.1.1:

    ```
    pip list | grep vos
    ```
    
    Upgrade vos if necessary:

    ```
    pip install --upgrade vos
    ```

    Uninstall `enum34`, the package raisng the AttributeError:
    
    ```
    pip uninstall enum34
    ```

    Then retry the `caom2utils` install:
    
    ```
    pip install caom2utils
    ```
