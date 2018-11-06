#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <thread>

#include "utility.h"

using namespace std;

void execute(string command) {
  int status = system(command.c_str());
  if (status < 0) {
    cerr << "Error: error while executing the command: " << command << endl;
    exit(-1);
  }
}

int main(int argc, char
  const * * argv) {

  if (argc != 4) {
    cerr << "USAGE: genesis INPUT.BAM REFERENCE.FA OUTPUT_DIRECTORY\n";
    exit(-1);
  }

  string input_bam = argv[1];
  string reference_genome = argv[2];
  string output_directory = argv[3];

  string input_bam_base_name_with_extension = input_bam.substr(input_bam.find_last_of("\\/") + 1);
  string str(input_bam_base_name_with_extension);

  string input_bam_base_name = str.substr(0, str.find_last_of("."));
  
/*
  //create output directory
  string create_output_directory = "mkdir -p " + output_directory;
  execute(create_output_directory);
*/

/*
  //split the input BAM file into mapped and unmapped bam file
  string split_to_mapped_bam_command = "../ext/samtools view -b -F 0x04 " + input_bam + " > " + output_directory + "/" + input_bam_base_name + "_mapped-reads.bam";
  string split_to_unmapped_bam_command = "../ext/samtools view -b -f 0x04 " + input_bam + " > " + output_directory + "/" + input_bam_base_name + "_unmapped-reads.bam";

  thread split_to_mapped_bam(execute, split_to_mapped_bam_command);
  thread split_to_unmapped_bam(execute, split_to_unmapped_bam_command);

  split_to_mapped_bam.join(); // pauses until split_to_mapped_bam finishes
  split_to_unmapped_bam.join(); // pauses until split_to_unmapped_bam finishes

  //generate statistics for the mapped and unmapped reads
  string generate_statistics_for_mapped_read_command = "../ext/samtools flagstat " + output_directory + "/" + input_bam_base_name + "_mapped-reads.bam > " + output_directory + "/" + input_bam_base_name + "_original_mapped_read.stat";
  string generate_statistics_for_unmapped_read_command = "../ext/samtools flagstat " + output_directory + "/" + input_bam_base_name + "_unmapped-reads.bam > " + output_directory + "/" + input_bam_base_name + "_original_unmapped_read.stat";

  thread generate_statistics_for_mapped_read(execute, generate_statistics_for_mapped_read_command);
  thread generate_statistics_for_unmapped_read(execute, generate_statistics_for_unmapped_read_command);

  generate_statistics_for_mapped_read.join(); // pauses until stat for mapped reads are generated
  generate_statistics_for_unmapped_read.join(); // pauses until stat for unmapped reads are generated

  //get reads from the originally unmapped bam file
  string get_reads_from_originally_unmapped_bam = "../ext/samtools bam2fq " + output_directory + "/" + input_bam_base_name + "_unmapped-reads.bam > " + output_directory + "/" + input_bam_base_name + "_unmapped-reads.fq";
  execute(get_reads_from_originally_unmapped_bam);
  
  // write code to do the quality control of the unmapped read using trimmomatic
*/
  //map the trimmed reads from the originally unmapped bam file to the reference using BWA-MEM
  string map_originally_unampped_read_to_reference_genome = "../ext/bwa mem -t 32 " + reference_genome + " " + output_directory + "/" + input_bam_base_name + "_unmapped-reads-trimmed-TruSeq2-SE.fq | ../ext/samtools view -bT " + reference_genome + " - > " + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.bam";
  execute(map_originally_unampped_read_to_reference_genome);

  //generate statistics about the newly mapped reads which were originally unmapped
  string generate_statistics_for_originally_unmapped_newly_mapped_reads = "../ext/samtools flagstat " + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.bam > " + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.stat";
  execute(generate_statistics_for_originally_unmapped_newly_mapped_reads);

  //sort the newly mapped bam file
  string sort_newly_mapped_bam = "../ext/samtools sort " + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.bam -o " + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted.bam";
  execute(sort_newly_mapped_bam);

  //mark PCR duplicates of the newly mapped bam file
  string mark_pcr_duplicate = "java -jar ../ext/MarkDuplicates.jar INPUT=" + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted.bam" + " OUTPUT=" + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup.bam" + " METRICS_FILE=" + output_directory + "/" + input_bam_base_name + "_metrics.txt";
  execute(mark_pcr_duplicate);

  //index the originally unmapped newly mapped sorted dedup bam file
  string index_merged_bam_file = "../ext/samtools index " + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup.bam";
  execute(index_merged_bam_file);

  //call snp and indel from the originally unmapped newly mapped sorted dedup bam file using Platypus
  string platypus_command = "python ../ext/Platypus_0.7.9.1/Platypus.py callVariants --bamFiles=" + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup.bam" + " --refFile=" + reference_genome + " --output=" + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup_variants.vcf";
  execute(platypus_command);
/*
  //extract indels from the Platypus generated vcf files
  generate_vcf_with_only_indel(output_directory + "/" + input_bam_base_name + "_merged_originally_mapped_overlapped_with_gene_region_and_newly_mapped_variants.vcf");

  //keep the genes that overlap with the indels called by Platypus
  string indel_set = input_bam_base_name + "_merged_originally_mapped_overlapped_with_gene_region_and_newly_mapped_variants_indel.vcf";
  check_overlap_between_gene_region_and_indel_set(gene_region_base_name, indel_set, input_bam_base_name, output_directory);
*/
  /*
    to be done later.
    -Keep only the indels that are in the gene_region / Keep only the genes that overlap with the indels
    -Delete all the temp files. Only keep the stat files and the snp and indel files.
  */
}
