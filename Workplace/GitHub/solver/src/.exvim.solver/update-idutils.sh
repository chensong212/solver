#!/bin/bash
export DEST="./.exvim.solver"
export TOOLS="/home/song/exvim/main//vimfiles/tools/"
export TMP="${DEST}/_ID"
export TARGET="${DEST}/ID"
sh ${TOOLS}/shell/bash/update-idutils.sh
