#! /usr/bin/bash
#
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]
then
	echo "Script is being sourced"
else
	echo "ERROR: Script is a subshell"
	echo "To affect your current shell enviroment source this script with:"
	echo "source env.sh"
	exit
fi

# check that this file is run
export ENV_RUN=true

# if acados folder not specified assume parent of the folder of the single examples
BLASFEO_INSTALL_DIR=${BLASFEO_INSTALL_DIR:-"$(pwd)/../blasfeo"}
export ACADOS_INSTALL_DIR
echo
echo "BLASFEO_INSTALL_DIR=$BLASFEO_INSTALL_DIR"

# if model folder not specified assume this folder
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BLASFEO_INSTALL_DIR/lib
echo
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
