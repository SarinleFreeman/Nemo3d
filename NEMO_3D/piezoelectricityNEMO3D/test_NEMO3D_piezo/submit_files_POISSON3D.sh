#!/bin/sh

# Submits all the named files as nemo3d PBS jobs

if [ x${1}x == "-h" ]; then
	cat <<EOF
Usage: submit_things.sh <template> <input1> <input2> <input3> <input4>
EOF
	exit 0
fi

template=$1
shift
file1=$1
shift
file2=$1
shift
file3=$1
shift
file4=$1

echo RUNDIR=potential_2_ham.d
echo "creating directory for the potential_2_ham"
echo if [ ! -d $RUNDIR ] then
echo	mkdir $RUNDIR
echo fi
mkdir -p potential_2_ham.d/
cp $file1 potential_2_ham.d/
cp $file2 potential_2_ham.d/
cp $file3 potential_2_ham.d/
cp $file4 potential_2_ham.d/
echo "Submitting $file1 $file2 $file3 $file4"
cat $template | sed -e "s/MY_INPUT_FILE_1/$file1/g" | sed -e "s/MY_INPUT_FILE_2/$file2/g" | sed -e "s/MY_INPUT_FILE_3/$file3/g" | sed -e "s/MY_INPUT_FILE_4/$file4/g" > ./potential_2_ham.d/potential_2_ham.pbs
# Remove the word 'echo' to make this go
cd potential_2_ham.d/
qsub potential_2_ham.pbs
cd ../..
