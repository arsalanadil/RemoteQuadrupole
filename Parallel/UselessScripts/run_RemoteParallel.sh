#!/bin/sh

max=2
for (( i=0; i < max; i++ ));
do
num=$((i + 1))
ssh "physics$num"
bash
cd ~/Parallel
python sqrt.py $num
done
