#!/bin/bash

docker run -it --rm --name caom2utilsTest caom2utils-363-test python caom2/setup.py test

docker run -v /home/goliaths/work/cadc/caom2tools/caom2utils:/usr/src/app/caom2utils --rm --name caom2utilsTest caom2utils-363-test
