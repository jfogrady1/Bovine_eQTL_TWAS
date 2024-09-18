#!/usr/bin/env python3

# This script will clean up a featureCounts matrix ouput
# Specifically, it will remove the long headings in sample names aswell as the .bam extension
# It will also remove the geneid name by modifying a temporary file which we will delete

import subprocess
import re
import sys
file = sys.argv[1]
tempfile = sys.argv[2]
outfile = sys.argv[3]


fin = open(file, 'r') # Open in read mode
ftemp = open(tempfile, 'w')

    

##########################
##Changing the headers####
##########################

# Extract the first line in input file
# Replace the file path with nothing in patterns
# Replace the .bam extension
# Write to the specified temp file

header = fin.readline()
header = header.replace('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/', '')
header= header.replace('_Aligned.sortedByCoord.out.bam', '')
ftemp.write(header)



#####################################
##Filtering out non expressed reads##
#####################################

# Cycle through lines in the file
# If the sum of the reads is < 2x the sample size
# Do not write them to the new file

for line in fin :
    fields = line.strip().split("\t")
    reads = 0
    if reads >= 0 : # this is what is modified for filtration
        ftemp.write(line)

  
ftemp.close()
fin.close()

################################
##Formating the geneID column###
################################

ftemp = open(tempfile, 'r') # open in read mode so as to not to overwrite the file
fout = open(outfile, 'w') # Open in write mode to write to write to at the end

header = ftemp.readline()
fout.write(header)

# Cycle through lines in the new file
# Change the first field = field[0] to the desired format of just GeneID:xxxxxxx
# Join the fields back together
# write to the file

for line in ftemp:
    fields = line.strip().split("\t")
    ID = fields[0]
    match = re.search(r'GeneID:\d*|Geneid', ID, re.I) 
    if match:
      new_id = match.group(0)
      new_id = new_id.replace(",", "")
      fields[0] = new_id
      new_line = '\t'.join(fields)
      fout.write(new_line + "\n")



fin.close()
fout.close()


      
