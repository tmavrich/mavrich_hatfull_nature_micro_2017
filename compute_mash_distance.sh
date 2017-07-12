#! /bin/bash
#Bash script to compute pairwise genomic distances using Mash
#Travis Mavrich
#version 2
#20150708





#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
if ! [[ $1 && $2 ]]
    then
        echo ""
        echo "This is a bash script to compute Mash genomic distances"
        echo ""
        echo "Execute script in working directory containing 'input_data' folder of fasta files"
        echo ""
        echo "First argument: sketch size (integer)"
        echo "Second argument: kmer size (integer)"
        echo ""
        echo "Outputs the mash distances file."
        echo "It also outputs a file containing list of all phage names used in the analysis."
        echo ""
        echo ""
        echo ""
        echo ""

        exit
fi




#Create the output file
sketch=$1
kmer=$2
currentdate=`date +%Y%m%d`
output_file=${currentdate}_${sketch}sketch_${kmer}kmer_mash.txt
working_dir=`pwd`
echo $kmer

#Create directories to store final files

output_dir=${working_dir}/${currentdate}_${sketch}sketch_${kmer}kmer_mash
sketch_dir=${output_dir}/sketch_output
mkdir $output_dir
mkdir $sketch_dir


#First create the sketches of each genome
echo Computing sketches...
for fasta_file in `ls ${working_dir}/input_data/`
do
    genome_name="${fasta_file%.*}"
    echo $genome_name >> ${currentdate}_${sketch}sketch_${kmer}kmer_mash_names.csv
	mash sketch -s $sketch -k $kmer ${working_dir}/input_data/${fasta_file}
	mv ${working_dir}/input_data/${fasta_file}.msh $sketch_dir
	
done



#Next compute distances from pairwise comparisons of sketches
echo Computing distances...
for first_sketch in `ls $sketch_dir`
do
	for second_sketch in `ls $sketch_dir`
	do	
		mash dist ${sketch_dir}/${first_sketch} ${sketch_dir}/${second_sketch} >> $output_file
	done	
done

mv $output_file $output_dir
