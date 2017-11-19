#!/bin/bash
#
# Initialize the environment for pwa_omega
#

# make sure ROOTSYS is defined
if ! echo $ROOTSYS | grep -q "/usr/local/root"; then
    export ROOTSYS=/usr/local/root
fi

# put root into the path, if not already there
if ! echo $PATH | grep -q "$ROOTSYS/bin"; then
    export PATH=$PATH:$ROOTSYS/bin
fi

# put root into the library path, if not already there
if ! echo $LD_LIBRARY_PATH | grep -q "$ROOTSYS/lib"; then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`root-config --libdir`
fi

# put root into the python path, if not already there
if ! echo $PYTHONPATH | grep -q "$ROOTSYS/lib"; then
    export PYTHONPATH=$PYTHONPATH:`root-config --libdir`
fi
