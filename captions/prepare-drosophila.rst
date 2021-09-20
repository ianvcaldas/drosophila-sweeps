
Turn the raw DGRP2 VCF files into data that we will apply machine learning models to.

The DGRP2 data we use comes from the "VCF file for the DGRP 2.0 Freeze", located `here <http://dgrp2.gnets.ncsu.edu/data.html>`_.

First, we remove every SNP in the data with more than 15% missing data. Then we perform imputation using Beagle.

From data data, we create lists of the DGRP lines with at least one resistant allele for each of the three loci of interest, as well as heterozygosity and allele counts for every site in the imputed dataset. Finally, we perform a selection scan where we calculate a variety of summary statistics across SNP windows of chromosomes 2R and 3R, where our known adaptive loci are located.
