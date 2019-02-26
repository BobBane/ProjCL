#!/bin/bash -e
# Pastes the common ProjCL stuff onto a checked-in SPA to create a shippable SPA
#
# Takes one argument - the destination directory
SPADIR=$1

# First, find out where we are

case $0 in
         /*)  SHELLFILE=$0 ;;
        ./*)  SHELLFILE=${PWD}${0#.} ;;
        ../*) SHELLFILE=${PWD%/*}${0#..} ;;
          *)  SHELLFILE=$(type -P $0) ; if [ ${SHELLFILE:0:1} != "/" ]; then SHELLFILE=${PWD}/$SHELLFILE ; fi ;;
esac
SHELLDIR=${SHELLFILE%/*}

JARS=jni/projcl_jni.jar
SLIBS="jni/libprojcl_jni_c.so src/libprojcl.so"

for x in $JARS
do
    cp ${SHELLDIR}/$x $SPADIR/algorithm/h2g/lib64
done

for x in $SLIBS
do
    cp ${SHELLDIR}/$x $SPADIR/algorithm/h2g/lib64/proj
done
