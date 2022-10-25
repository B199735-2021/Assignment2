# BPSM Assignment 2: A Python3 workflow

Due date: 1400 GMT, Mon 22 November 2021  
## Tasks overview
Write a generic Python3 programme/script 
  - that doesn't use BioPython
  - that isn't just a Unix/bash script called from Python3...
  - that you will make available as a passworded zip file in a PUBLIC repository on GitHub
  - that will allow the user:
  1. to identify a family of protein sequences from a user-defined subset of the taxonomic tree (e.g. glucose-6-phosphatase proteins from Aves (birds), or ABC transporters in mammals, or kinases in rodents, or adenyl
cyclases in vertebrates etc.) that could then be processed using, for example, one or more of the EMBOSS programmes installed on the MSc server:  
      - to determine, and plot, the level of protein sequence conservation across the species within that taxonomic group
      - to scan the protein sequence(s) of interest with motifs from the PROSITE database, to determine whether any known motifs (domains) are associated with this subset of sequences
      - to do any other appropriate EMBOSS (or other) analysis that you think might add relevant biological information to the outputs
  2. write a "help manual" for your programme which has two sections:
      - one aimed at an "ordinary user" (non-coder)
      - one aimed at a competent Python3 code-writer
