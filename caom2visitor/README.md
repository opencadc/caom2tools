# DOCUMENTATION

The caom2visitor script is a tool for iterating through the observations in a collection in a CAOM2 repo and applying changes to each observation according to an algorithm supplied by the user in a plugin file.

At a minimum, the plugin file should look like as follows:

   from caom2.caom2_observation import Observation

   class ObservationUpdater:
    
    def update(self, observation):
        assert isinstance(observation, Observation), (
            'observation {} is not an Observation'.format(observation))
        # custom code to update the observation

## System Requirments

* python2.7 or later

## Installation


To install caom2visitor, retrieve the [github](github.com/opencadc/caom2tools) distribution and use

 `python setup.py install --user`




## Development
A virtual environment (**venv**) is recommended to set up external dependencies without installing them system-wide. Following [these instructions](http://docs.python-guide.org/en/latest/dev/virtualenvs/), install **virtualenv**:
```
$ pip install virtualenv
```

Next, create, and activate a local **venv** (this example uses **bash**):
```
$ virtualenv venv
$ source venv/bin/activate
```
Use **pip** to install external dependencies used for testing this project (Note: if the following install fails you probably have missing/older version of required modules. In that case just run '$ pip install --upgrade <module>' and re-invoke the install command):
```
$ pip install -r dev_requirements.txt
```

After successfully installing the external dependencies, the unit tests are invoked by running
```
$ py.test [-v]
```

To obtain a coverage report use the following commands:
```
$ py.test --cov-report term-missing --cov=caom2visitor [-v] #command line report which includes missing lines

or

py.test --cov-report html:cov_html --cov=caom2visitor [-v] # html reports in the cov_html directory
```

Each time you resume work on the project and want to use the **venv** (e.g., from a new shell), simply re-activate it:
```
$ source venv/bin/activate
```
When done, just issue a 
```
$ deactivate
```
command to deactivate the virtual environment.

