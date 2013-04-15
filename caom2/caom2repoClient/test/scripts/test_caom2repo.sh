#!/bin/bash
#/*+
#************************************************************************
#****  C A N A D I A N   A S T R O N O M Y   D A T A   C E N T R E  *****
#*
#* (c) 2012.                            (c) 2012.
#* National Research Council            Conseil national de recherches
#* Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#* All rights reserved                  Tous droits reserves
#*
#* NRC disclaims any warranties,        Le CNRC denie toute garantie
#* expressed, implied, or statu-        enoncee, implicite ou legale,
#* tory, of any kind with respect       de quelque nature que se soit,
#* to the software, including           concernant le logiciel, y com-
#* without limitation any war-          pris sans restriction toute
#* ranty of merchantability or          garantie de valeur marchande
#* fitness for a particular pur-        ou de pertinence pour un usage
#* pose.  NRC shall not be liable       particulier.  Le CNRC ne
#* in any event for any damages,        pourra en aucun cas etre tenu
#* whether direct or indirect,          responsable de tout dommage,
#* special or general, consequen-       direct ou indirect, particul-
#* tial or incidental, arising          ier ou general, accessoire ou
#* from the use of the software.        fortuit, resultant de l'utili-
#*                                      sation du logiciel.
#*
#************************************************************************
#*
#*   Script Name:                       caom2repo.py
#*
#****  C A N A D I A N   A S T R O N O M Y   D A T A   C E N T R E  *****
#************************************************************************
#-*/

PYTHON=$(/usr/bin/which python)
TEST_SUBJECT="../../caom2repo.py -v"
VERSION_OUTPUT_1=$(${TEST_SUBJECT} --version 2>&1)

if [[ ! $VERSION_OUTPUT_1 = 'caom2repo.py 1.0' ]]
then
  echo "Failed version test."
  echo "Expected 'caom2repo.py 1.0', but got '${VERSION_OUTPUT_1}'"
  exit -1
fi

HELP_OUTPUT_1=$(${TEST_SUBJECT} 2>&1 | head -3 | tail -1)

if [[ ! $HELP_OUTPUT_1 = 'usage: caom2repo.py [-h] [--version] [-v] [-d]' ]]
then
  echo "Failed help test (1)."
  echo "Expected 'usage: caom2repo.py [-h] [--version] [-v] [-d]', but got '${HELP_OUTPUT_1}'"
  exit -1
fi

HELP_OUTPUT_2=$(${TEST_SUBJECT} -h 2>&1 | head -1)

if [[ ! $HELP_OUTPUT_2 = 'usage: caom2repo.py [-h] [--version] [-v] [-d]' ]]
then
  echo "Failed help test (2)."
  echo "Expected 'usage: caom2repo.py [-h] [--version] [-v] [-d]', but got '${HELP_OUTPUT_2}'"
  exit -1
fi

HELP_OUTPUT_3=$(${TEST_SUBJECT} --help 2>&1 | head -1)

if [[ ! $HELP_OUTPUT_3 = 'usage: caom2repo.py [-h] [--version] [-v] [-d]' ]]
then
  echo "Failed help test (3)."
  echo "Expected 'usage: caom2repo.py [-h] [--version] [-v] [-d]', but got '${HELP_OUTPUT_3}'"
  exit -1
fi

# Get tests
#
GET_OUTPUT_1=$(${TEST_SUBJECT} -g 2>&1 | tail -1)

if [[ ! $GET_OUTPUT_1 = 'caom2repo.py: error: argument -g/--get: expected 2 argument(s)' ]]
then
  echo "Failed get test (1)."
  echo "Expected 'caom2repo.py: error: argument -g/--get: expected 2 argument(s)', but got '${GET_OUTPUT_1}'"
  exit -1
fi

GET_OUTPUT_2=$(${TEST_SUBJECT} --get 2>&1 | tail -1)

if [[ ! $GET_OUTPUT_2 = 'caom2repo.py: error: argument -g/--get: expected 2 argument(s)' ]]
then
  echo "Failed get test (2)."
  echo "Expected 'caom2repo.py: error: argument -g/--get: expected 2 argument(s)', but got '${GET_OUTPUT_2}'"
  exit -1
fi

GET_OUTPUT_3=$(${TEST_SUBJECT} -g caom:TEST/BOGUSOBSERVATION 2>&1 | tail -1 2>&1)

if [[ ! $GET_OUTPUT_3 = 'caom2repo.py: error: argument -g/--get: expected 2 argument(s)' ]]
then
  echo "Failed get test (3)."
  echo "Expected 'caom2repo.py: error: argument -g/--get: expected 2 argument(s)', but got '${GET_OUTPUT_3}'"
  exit -1
fi

${TEST_SUBJECT} -g caom:IRIS/f212h000 /tmp/out.xml

GET_OUTPUT_4=$(cat /tmp/out.xml | sed 's/\r//' | head -1)
GET_OUTPUT_4_2=$(cat /tmp/out.xml | sed 's/\r//' | head -5 | tail -1)

if [[ ! ${GET_OUTPUT_4} = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" ]]
then
  echo "Failed get test (4)."
  echo "Expected '<?xml version=\"1.0\" encoding=\"UTF-8\"?>', but got '${GET_OUTPUT_4}'"
  exit -1
else
  if [[ ! ${GET_OUTPUT_4_2} = '  <caom2:metaRelease>2004-12-06T00:00:00.000</caom2:metaRelease>' ]]
  then
    echo "Failed get test (4.2)."
    echo "Expected '  <caom2:metaRelease>2004-12-06T00:00:00.000</caom2:metaRelease>', but got '${GET_OUTPUT_4_2}'"
    exit -1
  fi
fi

# Quick remove before a PUT
echo ""
echo "Issuing clean up delete in case test failed before final deletion.  Ignore any error message here."
${TEST_SUBJECT} -r caom:TEST/TESTBOGUS_CAOM2REPOCLIENT
echo "OK, pay attention again."

# Create tests
#
PUT_OUTPUT_1=$(${TEST_SUBJECT} -p 2>&1 | tail -1)

if [[ ! $PUT_OUTPUT_1 = 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)' ]]
then
  echo "Failed put test (1)."
  echo "Expected caom2repo.py: error: argument -p/--put: expected 2 argument(s)', but got '${PUT_OUTPUT_1}'"
  exit -1
fi

PUT_OUTPUT_2=$(${TEST_SUBJECT} --put 2>&1 | tail -1)

if [[ ! $PUT_OUTPUT_2 = 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)' ]]
then
  echo "Failed put test (2)."
  echo "Expected 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)', but got '${PUT_OUTPUT_2}'"
  exit -1
fi

PUT_OUTPUT_3=$(${TEST_SUBJECT} -p caom:TEST/TESTOBS 2>&1 | tail -1)

if [[ ! $PUT_OUTPUT_3 = 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)' ]]
then
  echo "Failed put test (3)."
  echo "Expected 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)', but got '${PUT_OUTPUT_3}'"
  exit -1
fi

PUT_OUTPUT_4=$(${TEST_SUBJECT} --put caom:TEST/TESTOBS 2>&1 | tail -1)

if [[ ! $PUT_OUTPUT_4 = 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)' ]]
then
  echo "Failed put test (4)."
  echo "Expected 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)', but got '${PUT_OUTPUT_4}'"
  exit -1
fi

PUT_OUTPUT_5=$(${TEST_SUBJECT} -p ../data/simple.xml 2>&1 | tail -1)

if [[ ! $PUT_OUTPUT_5 = 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)' ]]
then
  echo "Failed put test (5)."
  echo "Expected 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)', but got '${PUT_OUTPUT_5}'"
  exit -1
fi

PUT_OUTPUT_6=$(${TEST_SUBJECT} --put ../data/simple.xml 2>&1 | tail -1)

if [[ ! $PUT_OUTPUT_6 = 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)' ]]
then
  echo "Failed put test (6)."
  echo "Expected 'caom2repo.py: error: argument -p/--put: expected 2 argument(s)', but got '${PUT_OUTPUT_6}'"
  exit -1
fi

PUT_OUTPUT_7=$(${TEST_SUBJECT} -p caom:TEST/TESTBOGUS_CAOM2REPOCLIENT ../data/simple.xml 2>&1 | head -4 | tail -1)

if [[ ! $PUT_OUTPUT_7 = 'INFO:root:Successfully created Observation' ]]
then
  echo "Failed put test (7)."
  echo "Expected 'INFO:root:Successfully created Observation', but got '${PUT_OUTPUT_7}'"
  exit -1
fi

# Update tests
#
POST_OUTPUT_1=$(${TEST_SUBJECT} -u 2>&1 | head -1)

if [[ ! $POST_OUTPUT_1 = 'usage: caom2repo.py [-h] [--version] [-v] [-d]' ]]
then
  echo "Failed post test (1)."
  echo "Expected 'usage: caom2repo.py [-h] [--version] [-v] [-d]', but got '${POST_OUTPUT_1}'"
  exit -1
fi

POST_OUTPUT_2=$(${TEST_SUBJECT} --update 2>&1 | head -1)

if [[ ! $POST_OUTPUT_2 = 'usage: caom2repo.py [-h] [--version] [-v] [-d]' ]]
then
  echo "Failed post test (2)."
  echo "Expected 'usage: caom2repo.py [-h] [--version] [-v] [-d]', but got '${POST_OUTPUT_2}'"
  exit -1
fi

POST_OUTPUT_3=$(${TEST_SUBJECT} -u caom:TEST/TESTBOGUS_CAOM2REPOCLIENT 2>&1 | tail -1)

if [[ ! $POST_OUTPUT_3 = 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)' ]]
then
  echo "Failed post test (3)."
  echo "Expected 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)', but got '${POST_OUTPUT_3}'"
  exit -1
fi

POST_OUTPUT_4=$(${TEST_SUBJECT} --update caom:TEST/TESTOBS 2>&1 | tail -1)

if [[ ! $POST_OUTPUT_4 = 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)' ]]
then
  echo "Failed post test (4)."
  echo "Expected 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)', but got '${POST_OUTPUT_4}'"
  exit -1
fi

POST_OUTPUT_5=$(${TEST_SUBJECT} -u ../data/simple.xml 2>&1 | tail -1)

if [[ ! $POST_OUTPUT_5 = 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)' ]]
then
  echo "Failed post test (5)."
  echo "Expected 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)', but got '${POST_OUTPUT_5}'"
  exit -1
fi

POST_OUTPUT_6=$(${TEST_SUBJECT} --update ../data/simple.xml 2>&1 | tail -1)

if [[ ! $POST_OUTPUT_6 = 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)' ]]
then
  echo "Failed post test (6)."
  echo "Expected 'caom2repo.py: error: argument -u/--update: expected 2 argument(s)', but got '${POST_OUTPUT_6}'"
  exit -1
fi

POST_OUTPUT_7=$(${TEST_SUBJECT} -u caom:TEST/TESTBOGUS_CAOM2REPOCLIENT ../data/simple.xml 2>&1  | head -4 | tail -1)

if [[ ! $POST_OUTPUT_7 = 'INFO:root:Successfully updated Observation' ]]
then
  echo "Failed post test (7)."
  echo "Expected 'INFO:root:Successfully updated Observation', but got '${POST_OUTPUT_7}'"
  exit -1
fi

# Delete tests.
#
REMOVE_OUTPUT_1=$(${TEST_SUBJECT} -r caom:TEST/NOSUCHOBS 2>&1 | head -3 | tail -1)

if [[ ! $REMOVE_OUTPUT_1 = "ERROR:root:No such Observation found with URI 'caom:TEST/NOSUCHOBS'." ]]
then
  echo "Failed remove test (1)."
  echo "Expected 'ERROR:root:No such Observation found with URI 'caom:TEST/NOSUCHOBS'.', but got '${REMOVE_OUTPUT_1}'"
  exit -1
fi

REMOVE_OUTPUT_2=$(${TEST_SUBJECT} --remove caom:TEST/NOSUCHOBS 2>&1 | head -3 | tail -1)

if [[ ! $REMOVE_OUTPUT_2 = "ERROR:root:No such Observation found with URI 'caom:TEST/NOSUCHOBS'." ]]
then
  echo "Failed remove test (2)."
  echo "Expected 'ERROR:root:No such Observation found with URI 'caom:TEST/NOSUCHOBS'.', but got '${REMOVE_OUTPUT_2}'"
  exit -1
fi

REMOVE_OUTPUT_3=$(${TEST_SUBJECT} --remove caom:TEST/TESTBOGUS_CAOM2REPOCLIENT 2>&1 | head -3 | tail -1)

if [[ ! $REMOVE_OUTPUT_3 = "INFO:root:Successfully removed Observation caom:TEST/TESTBOGUS_CAOM2REPOCLIENT" ]]
then
  echo "Failed remove test (3)."
  echo "Expected 'INFO:root:Successfully removed Observation caom:TEST/TESTBOGUS_CAOM2REPOCLIEN', but got '${REMOVE_OUTPUT_3}'"
  exit -1
fi

echo "All tests passed."
