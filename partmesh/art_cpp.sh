#! /bin/bash

for i in *
do
# echo "Item $((l++)) : $i"

if [ -d $i ]; then
  cd $i
    echo "$PWD"
    astyle.exe --style=allman -s3 -C -S -K  *.cpp 
    astyle.exe --style=allman -s3 -C -S -K  *.h
  cd ..
fi

done
