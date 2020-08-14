#!/bin/bash
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License (GPL) version 3
# Earliest editions: 24/4/2017; the latest edition: 13/6/2018

##### Script information ##########
display_usage(){
    echo "
    Run hetSNP_depthPlot.R in parallel on Helix or other SLURM-supported Linux systems.
    Expect input files to be named in the format: [sample name][suffix (including the filename extension)]
    
    Options (11):
        --sampleName: a text file for a list of sample names
        --hetSNPdir: directory of parsed VCF files of heterozygous SNPs. No end forward slash.
        --homSNPdir: directory of parsed VCF files of homozygous SNPs. No end forward slash.
        --hetSNPsuffix: suffix of VCF files of heterozygous SNPs. For example, '__hetSNPs_info.tsv' for 'NJST258__hetSNPs_info.tsv'.
        --homSNPsuffix: suffix of VCF files of homozygous SNPs. For example, '__homSNPs_info.tsv' for 'NJST258__homSNPs_info.tsv'.
        --delim: deliminator of input SNP tables. For example, '\t' or ','.
        --outdir: directory for output figures
        --figureSuffix: suffix inserted between the sample name and the filename extension (.png) of output figures
        --genomeLen: length of the reference genome
        --codeDir: directory of hetSNP_depthPlot.R without end forward slash
        --bundleSize: number of jobs per bundle
        
    Example:
        bash hetSNP_depthPlot_slurm.sh --sampleName sample_list.txt --hetSNPdir het --homSNPdir hom --hetSNPsuffix '__hetSNPs_info.tsv' \
        --homSNPsuffix '__homSNPs_info.tsv' --delim '\t' --outdir pic --figureSuffix hetSNPs --genomeLen 5.25e6 --codeDir codes \
        --bundleSize 16
        
    Output file: [outdir]/[number of hetSNPs]__[sample name]__[figure suffix].png
    "
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
    display_usage
    exit
fi

# constants and defaults
R_MODULE="R/3.3.3-vlsci_intel-2015.08.25"  # for Helix users
bundle_size=16
code_dir="/vlsci/SG0006/shared/wan/scripts/contamination_assessment"
PAR="sysgen"  # job partition

##### Main program ##########
# Read arguments
while [[ $# -gt 1 ]]; do  # loops as long as the number of arguments >= 2
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
        --figureSuffix)
        fig_suf="$2"
        shift 2
        ;;
        --genomeLen)
        genome_len="$2"
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
        --walltime)
        walltime="$2"
        shift 2
        ;;
        *)
        shift 2  # deliberately left blank for undefined arguments
        ;;
    esac  # end of the case statement
done

# Check whether the output directory exists or not
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

# Submit jobs
n_jobs=`cat ${samples} | wc -l`  # number of jobs
echo "Processing ${n_jobs} samples."

bundle_index=0
n=0
cmd=""
job_processed=0
walltime="0-1:00:00"

while IFS= read -r id ; do  # read every line of a text file
    # Skip strains that do not have any hetSNPs, which do not have the corresponding TSV files produced as well by filter_SNP_info.R.
    hetSNP_tab="${het_dir}/${id}${het_suf}"  # coordinate table for hetSNPs
    homSNP_tab="${hom_dir}/${id}${hom_suf}"
    if [ ! -f "$hetSNP_tab" ]; then
        printf "Skip ${id} because its hetSNP table ${hetSNP_tab} is absent.\n"
        job_processed=$(( job_processed+1 ))
        continue  # start the next iteration
    fi
    if [ ! -f "$homSNP_tab" ]; then
        printf "Skip ${id} because its homSNP table ${homSNP_tab} is absent.\n"
        job_processed=$(( job_processed+1 ))
        continue
    fi
    
    # Else, create a job and submit it to the SLURM system.
    cmd+="srun --nodes=1 --ntasks=1 --cpus-per-task=1 Rscript ${code_dir}/hetSNP_depthPlot.R --hetSNP ${hetSNP_tab} --homSNP ${homSNP_tab} --delim '${delim}' --sampleName ${id} --genomeLen ${genome_len} --outdir ${outdir} --suffix ${fig_suf} &\n"
    n=$(( n+1 ))
    job_processed=$(( job_processed+1 ))
    if [ "$n" -eq "$bundle_size" ] || [ "${job_processed}" -eq "${n_jobs}" ]; then  # The second statement accounts for the last bundle of jobs, which may be less than 16.
        bundle_index=$(( bundle_index + 1 ))  # wrap jobs into a new bundle
        script="bundle__${bundle_index}.slurm"  # bundle single jobs into a single script
        printf "#!/bin/bash\n" > $script  # The echo command does not interpret the newline character.
        printf "#SBATCH -p ${PAR}\n" >> $script
        printf "#SBATCH --ntasks=${bundle_size}\n" >> $script
        printf "#SBATCH --mem-per-cpu=512\n" >> $script
        printf "#SBATCH --time=${walltime}\n" >> $script
        printf "#SBATCH --job-name=bundle_${bundle_index}\n\n" >> $script
        printf "module load ${R_MODULE}\n\n" >> $script
        printf "cd ${PWD}\n\n" >> $script
        printf "${cmd}\n" >> $script
        printf "wait\n" >> $script  # IMPORTANT: must wait for all subordinate jobs to finish before finishing the parental job (namely, the bundle), or all get killed once the parental script finishes.
        sbatch $script
        n=0
        cmd=""
		sleep 1  # to allow the scheduler to process all commands
    fi
done < $samples

echo "Done!"
