#!/bin/bash
export DEST="./.exvim.solver"
export TOOLS="/home/song/exvim/main//vimfiles/tools/"
export TMP="${DEST}/_inherits"
export TARGET="${DEST}/inherits"
sh ${TOOLS}/shell/bash/update-inherits.sh
