#!/bin/bash
#
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License (GPL) version 3
# First edition: 22/4/2017; the latest edition: 12/6/2018

display_usage() {
    echo "
    Filter extracted information of VCFs in accordance with genomic coordinates specified in a text file.
    Useage:
      ./filterSNPs_slurm.sh -s [sample list] -i [info_dir] -o [out_dir] -f [filter file]
    Example:
      ./filterSNPs_slurm.sh -s sequence_list.txt -i cgSNPs/hom -o cgSNPs/hom/filtered -f Kp_rm.coords
    Notice directory names should not have the ending forward slash. In addition, filter_SNP_info.R and this
    script must be put in the same directory.
    "
}

# Display usage information if a null argument, "-h" or "--help" is encountered.
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	display_usage
	exit
fi

# Environmental configurations
PAR="sysgen"  # queue/paritition name
R_MODULE="R/3.3.3-vlsci_intel-2015.08.25"
script_dir=`dirname $0`  # get where the scripts are stored
wait_time=1

# Read arguments
while [[ $# -gt 1 ]]; do  # loops as long as the number of arguments >= 2
    key="$1"
    case $key in
        -i|--info_dir)
        info_dir="$2"
        shift 2 # remove the first two arguments from the vector; otherwise, the script enters a dead loop.
        ;;
        -s|--samples)
        samples="$2"
        shift 2
        ;;
        -o|--outdir)
        out_dir="$2"
        shift 2
        ;;
        -f|--filter)
        filter="$2"
        shift 2
        ;;
		-w|--wait)
		wait_time=$2  # time interval before job submissions
		shift 2
		;;
        *)
        shift 2  # deliberately left blank for undefined arguments
        ;;
    esac  # end of the case statement
done

# Check and set up the output directory
if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

# Submit a job for each sample
printf "\nMaking an SLURM script for every VCF file.\n"
while read -r id; do
	# skip the current sample if its output files are present
	if [ -e "${out_dir}/${id}__info_retained.tsv" ] && [ -e "${out_dir}/${id}__info_snpNum.tsv" ] && [ -e "${out_dir}/${id}__info_removed.tsv" ]; then
		echo "Skip the sample ${id} as all of its three output files are found."
	else
		script="${id}.slurm"
		printf "#!/bin/bash\n" > $script  # The echo command does not interpret the newline character.
		printf "#SBATCH -p ${PAR}\n" >> $script
		printf "#SBATCH --job-name='Filter_${id}'\n\n" >> $script
		printf "#SBATCH --ntasks=1\n" >> $script
		printf "#SBATCH --nodes=1\n" >> $script
		printf "#SBATCH --ntasks-per-node=1\n" >> $script
		printf "#SBATCH --cpus-per-task=1\n" >> $script
		printf "#SBATCH --mem-per-cpu=2048\n" >> $script
		printf "#SBATCH --time=0-5:0:00\n\n" >> $script
		printf "module load ${R_MODULE}\n\n" >> $script
		printf "Rscript ${script_dir}/filter_SNP_info.R ${info_dir}/${id}__info.tsv ${filter} ${out_dir}/${id}__info ${id}\n" >> $script
		sbatch $script
		sleep $wait_time  # to allow the scheduler to process all commands
	fi
done < "$samples"

echo "Finish!"