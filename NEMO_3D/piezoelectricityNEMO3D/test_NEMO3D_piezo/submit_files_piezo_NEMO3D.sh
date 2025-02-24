#!/bin/sh

# Submits all the named files as nemo3d PBS jobs

if [ x${1}x == "-h" ]; then
	cat <<EOF
Usage: submit_things.sh <template> <#nodes> <input1> <input2> ... <inputN>
EOF
	exit 0
fi

template=$1
shift
nodes=$1
shift
files=$1

for file in $files; do
	echo RUNDIR=piezo_files.d
	echo "creating directory for the piezo_files"
	echo if [ ! -d $RUNDIR ] then
	echo	mkdir $RUNDIR
	echo fi
	mkdir -p piezo_files.d/$nodes.d
        cp $file piezo_files.d/$nodes.d/
        cp potential piezo_files.d/$nodes.d/
	echo "Submitting $file"
	cat $template | sed -e "s/MY_NODES/$nodes/g" | sed -e "s/MY_INPUT_FILE/$file/g" > ./piezo_files.d/$nodes.d/piezo_files.$nodes.pbs
	# Remove the word 'echo' to make this go
	cd piezo_files.d/$nodes.d
	qsub piezo_files.$nodes.pbs
	cd ../..
done
