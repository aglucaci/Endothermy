#!/bin/bash

FASTA=$1
OUTPUT=$2

echo "# Input FASTA: "$FASTA
echo "# Output FASTA: "$OUTPUT

# Creating a copy of the data in the results folder
echo ""
echo "Creating a copy of the data, placed in the results folder"
echo "cp $FASTA ${OUTPUT}"
cp $FASTA ${OUTPUT}

#echo ""
#echo "Cleaning our FASTA header information in: "${OUTPUT}.fa
#echo "Replacing illegal characters with underscores"
#sed -i 's/ /_/g'  ${OUTPUT}
#sed -i 's/:/_/g'  ${OUTPUT}
#sed -i 's/|/_/g'  ${OUTPUT}
#sed -i 's/-/_/g'  ${OUTPUT}
#sed -i 's/\//_/g' ${OUTPUT}

exit 0