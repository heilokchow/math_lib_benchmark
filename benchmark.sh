#!/bin/bash

if [ -f "result.txt" ]
then
	rm result.txt
fi

if [ -f "a.out" ]
then 
	make clean
fi

for x in 1 2 3 4 5
do
	for y in 1000 2000 4000 8000
	do
		make MODEL=$x
		./a.out $y 5 >> log 2>&1
		make clean
	done
done		

date >> log 2>&1
