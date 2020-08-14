# Code for Detecting DNA Contamination in Sequence Reads

Yu Wan

First release: 13 Aug 2020; last modification: 14 Aug 2020



This repository offers code for visual evaluation of DNA contamination in sequence reads. Most code works on the variant call format (VCF). The latest specification is [VCF v4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf).



Dependencies:

- [VCFtools](https://vcftools.github.io/index.html)
- R, Rscript, and package optparse
- Linux bash



## 1. extractInfoFromVCF.sh

This bash script extracts values of sub-fields from field INFO using VCFtools and produces a tab-delimited file for each sample. Assuming every name of VCF files follows format `[sample name][shared root of filenames].vcf` (For example, `DT104x_hetSNPs.vcf`, where `DT104x` is the sample name and `_hetSNPs` is the filename root), this script can be run with the following parameters:

- `-i/--input`: A text file listing sample names for which VCF files will be accessed. One line per name.
- `-s/--source`: Directory where VCF files are stored.
- `-o/--output`: Output directory.
- `-r/--root`: Root shared by filenames of VCF files.
- `-f/--fields`: Names of sub-fields within field INFO for value extraction. Since field names POS, REF, ALT are not part of INFO and will always be extracted by VCFtools, they must not be placed in this parameter.
- `-d/--vcftoolsDir`: Installation directory of VCFtools. Default: current working directory.

Example command line:

`bash ~/cgSNPs/contamination/extractInfoFromVCF.sh -i sample_id.txt -s $PWD -o $PWD/het -r 'study1_sorted_qual_hetSNPs' -f 'DP,DP4' -d "${HOME}/anaconda3/envs/mapping/bin"`

Output filenames follow the format `[sample name]__info.tsv`. An example output (`DT104x__info.tsv`):

```text
CHROM	POS	REF	ALT	DP	DP4
DT104x_chr	62354	T	C	73	10,0,27,2
DT104x_chr	101232	C	G	45	15,2,3,2
DT104x_chr	103233	T	C	45	1,0,3,2
DT104x_chr	455942	G	T	47	2,2,22,2
```

Use command `which vcftools` to find out the installation directory of VCFtools.



## 2. hetSNP_depthPlot.R

(To be continue)



