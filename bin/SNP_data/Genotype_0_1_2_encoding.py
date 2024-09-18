# This script converts a vcf to 0/1/2 encoding as a csv format as required by Peer
# and Matrix eQTL


# Import relevant modules
import os
import subprocess
import re
import sys

file = sys.argv[1]
outfile = sys.argv[2]

# Open necessary files
genotype = open(file, "r")
genotype_peer = open(outfile, "w")

## Genotype first


# Cycle through the file
# If we meet the default header, split at tab
# in the new file write the ID
for line in genotype:
  if line.startswith("#CHROM"):
    fields = line.strip().split("\t")
    genotype_peer.write("ID" + "\t")
    
    
    # From field 9 to the end of the file (len of fields)
    # Replace everything if names have not been changed to simple names
    # This may have been done already
    for i in range(9, len(fields)):
      sample = fields[i]
      sample = sample.replace('"', "")
      sample = re.sub(r'[A-Z]\d_', '', sample)
      sample = re.sub(r'_.*', '', sample)
      genotype_peer.write(sample + "\t")
      if i == len(fields) - 1: # -1 as python is 0 based
        genotype_peer.write("\n")
  # For every other line, split at tab
  # recode the genotype 0/0 = 0, 0/1 = 1 , 1/1 = 2
  elif line.find("#") == -1:
    fields2 = line.strip().split("\t")
    ID = str(fields2[2])
    ID = ID.replace('"', "")
    genotype_peer.write(fields2[2] + "\t")
    for n in range(9, len(fields2)):
      geno = fields2[n]
      geno = geno.replace('"', "")
      geno = str(geno)
      if geno == "0|0":
        geno = "0"
        genotype_peer.write(geno + "\t")
      elif geno == "0|1":
        geno = "1"
        genotype_peer.write(geno + "\t")
      elif geno == "1|0":
        geno = "1"
        genotype_peer.write(geno + "\t")
      elif geno == "1|1":
        geno = "2"
        genotype_peer.write(geno + "\t")
      elif geno == ".|.":
        geno = "NA"
        genotype_peer.write(geno + "\t")
      if n == len(fields2) - 1:
        genotype_peer.write("\n")  # create a new line
        
# Close the file         
genotype_peer.close()


