#!/bin/bash

FILE=$1
OUTPUT=$2

#awk -F \"t\" '{if($15 == 0 && $16==1 && $5>=3 && $6>=0.95){print}}\' $1 | awk -F "\t" '{if($10>=30){print}}' | awk -F "\t" '{if($19!="NA"){print}}' > $OUTPUT

# Col 19 is panel of normal number of samples
# Change PON filter to 50 for merged cohort (initial == 1)
awk -F "\t" '{if($18 == 0 && $19<=50 && $5>=3 && $6>=0.95 && $10>=30 && $22!="NA"){print}}\' $FILE

#awk -F "\t" '{print $1}' $FILE


