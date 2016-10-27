import sys
import os
from caom2visitor import __version__

if sys.version_info[0] > 2:
    print "The  is only compatible with Python version 2.n"
    sys.exit(-1)

# Build the list of scripts to be installed.

dependencies = ['CAOM2RepoClient', 'lxml', 'pyCAOM2']

script_dir = 'scripts'
scripts = []
for script in os.listdir(script_dir):
    if script[-1] in ["~", "#"]:
        continue
    scripts.append(os.path.join(script_dir,script))

try:
    from setuptools import setup, find_packages
    has_setuptools = True
except:
    from distutils.core import setup
    has_setuptools = False


execfile('caom2visitor/__version__.py')

setup(name="caom2visitor",
      version=__version__.version,
      url="https://github.com/opencadc/caom2tools/caom2visitor",
      description="Tool to update the records in a CAOM2 data collection.",
      author="Adrian Damian",
      maintainer="CADC",
      maintainer_email="cadc@nrc-cnrc.gc.ca",
      long_description="""A tool that visits observations in a data collection.
                          It loads every visited observation and calls a user
                          provided python method to update it before
                          persisting it back in the collection.""",
      license="AGPLv3",
      packages=find_packages(exclude=['']),
      scripts=scripts,
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU Affero General Public License v3',
          'Operating System :: POSIX',
          'Programming Language :: Python',
      ],
      install_requires=dependencies
      )




