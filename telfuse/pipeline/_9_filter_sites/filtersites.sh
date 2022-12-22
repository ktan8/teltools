#!/bin/bash

FILE=$1
OUTPUT=$2

#awk -F \"t\" '{if($15 == 0 && $16==1 && $5>=3 && $6>=0.95){print}}\' $1 | awk -F "\t" '{if($10>=30){print}}' | awk -F "\t" '{if($19!="NA"){print}}' > $OUTPUT


awk -F "\t" '{if($15 == 0 && $16==1 && $5>=3 && $6>=0.95 && $10>=30 && $19!="NA"){print}}\' $FILE

#awk -F "\t" '{print $1}' $FILE


