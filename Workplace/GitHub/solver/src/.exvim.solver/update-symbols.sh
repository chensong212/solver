#!/bin/bash
export DEST="./.exvim.solver"
export TOOLS="/home/song/exvim/main//vimfiles/tools/"
export TMP="${DEST}/_symbols"
export TARGET="${DEST}/symbols"
sh ${TOOLS}/shell/bash/update-symbols.sh
