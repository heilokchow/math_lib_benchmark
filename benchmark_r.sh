#!/bin/bash

if [ -f "result.txt" ]
then
	rm result.txt
fi

if [ -f "a.out" ]
then 
	make clean
fi

for x in 3
do
	for y in 50 100 400 800 1200 1600 2000
	do
		make MODEL=$x
		./a.out $y 5 >> log 2>&1
		make clean
	done
done		

date >> log 2>&1
