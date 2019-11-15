#!/bin/bash -l

mkdir bin
cd src
search_dir=''
src_files=''
for entry in *.cpp
do
  src_files="$src_files"" ""$entry"
done

#  sbatch "$entry"

command="g++ -O3 -std=c++11 -lpthread -std=c++0x -pthread -I.$src_files -o ../bin/MC-DC_Simulator"


eval "$command"