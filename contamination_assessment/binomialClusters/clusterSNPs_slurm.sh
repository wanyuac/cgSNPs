#!/bin/bash

##### Script information ##########
display_usage(){
    echo "
    Run clusterSNPs.R in parallel on Helix or other SLURM-supported Linux systems.
    Expect input files to be named in the format: [sample name][suffix (including the filename extension)]
    
    Options:
        --sampleName: a text file for a list of sample names
        --hetSNPdir: directory of parsed VCF files of heterozygous SNPs. No end forward slash.
        --homSNPdir: directory of parsed VCF files of homozygous SNPs. No end forward slash.
        --hetSNPsuffix: suffix of VCF files of heterozygous SNPs. For example, '__hetSNPs_info.tsv' for 'NJST258__hetSNPs_info.tsv'.
        --homSNPsuffix: suffix of VCF files of homozygous SNPs. For example, '__homSNPs_info.tsv' for 'NJST258__homSNPs_info.tsv'.
        --delim: deliminator of input SNP tables. For example, '\t' or ','.
        --outdir: directory for the output file
        --outputSuffix: filename suffix of the final output
        --codeDir: directory of hetSNP_depthPlot.R without end forward slash
        --bundleSize: number of tasks per core
        --override: flag it to not skip jobs whose outputs have been generated
        --debug: flag it to avoid submitting jobs to SLURM. In this case, SLURM scripts will still be generated.
        
    Example:
        bash clusterSNPs_slurm.sh --sampleName strains_all.txt --hetSNPdir het --homSNPdir hom --hetSNPsuffix '__hetSNPs_info.tsv' \
        --homSNPsuffix '__homSNPs_info.tsv' --delim '\t' --outdir binommix --outputSuffix SNPclusters --codeDir hetSNP_inspector \
        --bundleSize 10 --debug --override # 
        
    Output file for each sample: [outdir]/[sample name]__[output suffix].tsv (a tab-delimited text file)
        
    Copyright 2017 Yu Wan <wanyuac@gmail.com>
    Licensed under the Apache License, Version 2.0
    Editions: 4 May 2017
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

# constants
R_MODULE="R/3.3.3-vlsci_intel-2015.08.25"  # for Helix users
PARTITION="sysgen"  # Which partition of job queues the SLURM script will be submitted to.
MEM="512"  # memory assigned to every core
WALLTIME_TOTAL="2-12:00:00"
WALLTIME_PERJOB="6:0:0"

# default argument
debug=false
skip=true

##### Main program ##########
# read arguments
while [[ $# -gt 0 ]]; do  # loops as long as the number of arguments >= 2
    key="$1"  # read the first argument in the argument string
    case $key in
        --sampleName)
        samples="$2"
        shift 2 # remove the first two arguments from the vector
        ;;
        --hetSNPdir)
        het_dir="$2"
        shift 2
        ;;
        --homSNPdir)
        hom_dir="$2"
        shift 2
        ;;
        --hetSNPsuffix)
        het_suf="$2"
        shift 2
        ;;
        --homSNPsuffix)
        hom_suf="$2"
        shift 2
        ;;
        --delim)
        delim="$2"
        shift 2
        ;;
        --outdir)
        outdir="$2"
        shift 2
        ;;
        --outputSuffix)
        out_suf="$2"
        shift 2
        ;;
        --codeDir)
        code_dir="$2"
        shift 2
        ;;
        --bundleSize)
        bundle_size="$2"
        shift 2
        ;;
        --override)
        skip=false
        shift
        ;;
        --debug)
        debug=true
        shift  # This option only consumes a single argument in the list.
        ;;
        *)
        shift  # deliberately left blank for undefined arguments; consume a single argument each time
        ;;
    esac  # end of the case statement
done

# make a SLURM script and submit jobs (when it is configured to do so)
sample_n=`cat ${samples} | wc -l`  # number of samples to be processed
sample_count=0  # number of samples processed
n_jobs=0
bundle_index=0
cmd=""

echo "Processing data of ${sample_n} samples."
while IFS= read -r id ; do  # read every sample name in the text file
    out_file="${outdir}/${id}__${out_suf}.tsv"
	sample_count=$(( sample_count+1 ))
    if $skip; then
        if [[ -e "${out_file}" ]]; then
            echo "Existing file ${out_file} is skipped."
        else
            cmd+="srun --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=${MEM} --time=${WALLTIME_PERJOB} Rscript ${code_dir}/clusterSNPs.R --hetSNP ${het_dir}/${id}${het_suf} --homSNP ${hom_dir}/${id}${hom_suf} --delim '${delim}' --sampleName ${id} --outdir ${outdir} --suffix ${out_suf}\n"
            n_jobs=$(( n_jobs+1 ))
        fi
    else
        cmd+="srun --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=${MEM} --time=${WALLTIME_PERJOB} Rscript ${code_dir}/clusterSNPs.R --hetSNP ${het_dir}/${id}${het_suf} --homSNP ${hom_dir}/${id}${hom_suf} --delim '${delim}' --sampleName ${id} --outdir ${outdir} --suffix ${out_suf}\n"
        n_jobs=$(( n_jobs+1 ))
    fi
	
	# make and submit a SLURM script once the number of jobs reaches the bundle size or the number of samples (for the last bundle)
	if [ "$n_jobs" -eq "$bundle_size" ] || [ "${sample_count}" -eq "${sample_n}" ]; then
		echo "Submitting ${n_jobs} jobs."
		bundle_index=$(( bundle_index+1 ))
		script="clusterSNPs_${bundle_index}.slurm"  # bundle single jobs into a single script
		printf "#!/bin/bash\n" > $script  # The echo command does not interpret the newline character.
		printf "#SBATCH --job-name=clusterSNPs_${bundle_index}\n\n" >> $script
		printf "#SBATCH -p ${PARTITION}\n" >> $script
		printf "#SBATCH --time=${WALLTIME_TOTAL}\n" >> $script
		printf "module load ${R_MODULE}\n\n" >> $script
		printf "cd ${PWD}\n\n" >> $script
		printf "${cmd}" >> $script
		if ! $debug; then
			sbatch $script  # allocate resources and submit a batch script
			sleep 1
		fi
		cmd=""
		n_jobs=0
	fi
done < $samples

echo "Done!"
