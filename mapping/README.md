# Pipeline for Read Mapping

Yu Wan

23 Aug 2020



This repository offers a basic bash-based pipeline for mapping short reads against a reference genome and call variants from BAM files. Optionally, it generates a depth plot about homozygous SNPs and heterozygous SNPs for each sample using R script `hetSNP_depthPlot.R` in subdirectory `contamination`. In order to run this pipeline, users need submit script `read_mapping.pbs` to the Portable Batch System (PBS) on a high-performance computing cluster (HPC). Currently, the pipeline only supports Bowtie 2.



Dependencies:

- Bowtie 2 (step 1)
- BCFtools (step 2, variant calling)
- VCFtools (step 3, plotting)
- Rscript (step 3, plotting)



The main script of this pipeline is `read_mapping.pbs`, which loads user configurations from text file `read_mapping.config`. Script `read_mapping.pbs` takes a single argument `c` for accessing the configuration file. An example command is shown as follows:

```bash
qsub -v c=read_mapping.config read_mapping.pbs
```

Please see configuration file `read_mapping.config` for parameters that need to be set by users before executing this pipeline.