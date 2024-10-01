#!/usr/bin/env bash

set -e

if [[ $# != 2 ]]
then
	echo "usage $0 <crs_matrix>"
else
	set -x

	MAT_PATH=$1
	EXEC=./partition_crs

	mpiexec -np 2 $EXEC $MAT_PATH/rowptr.raw $MAT_PATH/colidx.raw int int &

	sleep 1
	ps aux | grep $EXEC | grep rowptr | awk '{print $2}'| head -2 > process.txt
	bat process.txt

	proc_id=`head -1 process.txt`
	lldb -p $proc_id $EXEC 
fi
