# PLOS_One_e0194060_rebuttal

This repository was created by Garrett A. Meek (garrett.a.meek@gmail.com),
and contains files used to reproduce the analyses discussed in:

Carmody LA, Caverly LJ, Foster BK, Rogers MAM, Kalikin LM, et al. (2018) Fluctuations in airway bacterial communities associated with clinical states and disease stages in cystic fibrosis. PLOS ONE 13(3): e0194060.

In order to reproduce this protocol, simply execute 'run.csh' from the command line on any Linux OS:

./run.csh > run.out

The following external software is required to run this protocol:

SRA Toolkit (v 2.8.2-1)
EDirect (v 10.7)
mothur (v 1.39.5)

Please note that the provided version of 'run.csh' skips initial data gathering steps in its present form. This can be changed with:

set download_fastq_files = "True"
set get_new_data = "True"
