#!/bin/bash -e
# Creates an array of PLKernelSources elements, populated
# with filenames and string contents of files, hosed into
# projcl_kernel_sources_text.h in src.
# Won't work on Windows, sorry...

# Need the directory this script is in
case $0 in
         /*)  SHELLFILE=$0 ;;
        ./*)  SHELLFILE=${PWD}${0#.} ;;
        ../*) SHELLFILE=${PWD%/*}${0#..} ;;
          *)  SHELLFILE=$(type -P $0) ; if [ ${SHELLFILE:0:1} != "/" ]; then SHELLFILE=${PWD}/$SHELLFILE ; fi ;;
esac
SHELLDIR=${SHELLFILE%/*}

# So, SHELLDIR is projectHome/kernel.  We want to run from there...
cd $SHELLDIR

# and put our output into...
KSTFILE=../src/projcl_kernel_sources_text.h

# Clean out the old one
rm -f $KSTFILE

# Hose the data files into C string arrays
# Note the ./ in the loop - this insures the file names come out as
#     ./mumble.opencl
# which xxd will render in variable names as
# __mumble_opencl
# which shouldn't collide with anything else
for x in ./*.opencl
do
    xxd -i $x | sed 's/\([0-9a-f]\)$/\0, 0x00/' | sed 's/^unsigned int.*$//' \
	>> $KSTFILE
done

echo 'static PLKernelSources _pl_kernelsources[] = {' >> $KSTFILE

comma=
for x in ./*.opencl
do
    echo $comma '{ "'${x##*/}'"' ',' ${x//[\/.]/_} '}' >> $KSTFILE
    comma=,
done
# And have a last entry we can test for
echo ', { NULL, NULL }' >> $KSTFILE
echo '};' >> $KSTFILE
