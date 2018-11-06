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

void generate_vcf_with_only_indel(string rawVCFFile) {
  string rawVCFFileNameWithoutExtension = rawVCFFile.substr(0, rawVCFFile.find_last_of("."));
  string indel_file_name = rawVCFFileNameWithoutExtension + "_indel.vcf";

  ifstream fin(rawVCFFile);
  FILE * fout = fopen((indel_file_name).c_str(), "w");

  if (fout != NULL && fin.is_open()) {
    string line;
    while (getline(fin, line)) {
      if (line.at(0) == '#') {
        fprintf(fout, "%s\n", line.c_str()); // only the header of the VCF file
      } else {
        vector <string> tokens = split(line, '\t');
        std::size_t found;
        string chrom = tokens.at(0); // CHROM
        string pos = tokens.at(1); // POS
        string id = tokens.at(2); // ID
        string refer = tokens.at(3); // REF
        string alt = " "; // ALT will be filled later
        string qual = tokens.at(5); // QUAL
        string filter = tokens.at(6); // FILTER
        string info = tokens.at(7); // INFO
        string format = tokens.at(8); // FORMAT
        string sampleId = tokens.at(9); // SAMPLE_ID

        found = tokens.at(4).find(','); // some rows have more than one indel in the alt column
        if (found != string::npos) {
          vector <string> alternative_indel_token = split(tokens.at(4).c_str(), ',');
          for (int i = 0; i < alternative_indel_token.size(); i++) {
            if ((alternative_indel_token.at(i)).at(0) != '<') { // get rid of contents like <SV>
              alt = alternative_indel_token.at(i);
              if (refer.length() != alt.length()) {
                fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                  chrom.c_str(), pos.c_str(), id.c_str(), refer.c_str(),
                  alt.c_str(), qual.c_str(), filter.c_str(), info.c_str(), format.c_str(), sampleId.c_str());
              }
            }
          }
        } else {
          if (refer.length() != tokens.at(4).length()) {
            string alt = tokens.at(4);
            //fprintf(fout, "%s\n", line.c_str());
            fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
              chrom.c_str(), pos.c_str(), id.c_str(), refer.c_str(),
              alt.c_str(), qual.c_str(), filter.c_str(), info.c_str(), format.c_str(), sampleId.c_str());
          }
        }
      }
    }
    fclose(fout);
  } else {
    cerr << "ERROR: File operation error while extracting indels.\n";
    exit(-1);
  }
}
void check_overlap_between_gene_region_and_indel_set(string gene_region_base_name, string indel_set,
  string input_bam_base_name, string output_directory) {

  //check overlap between sorted gene_region and sorted indel_set
  string overlap_between_gene_region_and_indel_set = "../ext/bedtools_2.17.0 intersect -wa -wb -a " 
	+ output_directory + "/" + gene_region_base_name + "_sorted.bed" + " -b " 
	+ output_directory + "/" + indel_set + " > " 
	+ output_directory + "/" + input_bam_base_name + "_gene_regions_overlapped_with_indels.txt";
  execute(overlap_between_gene_region_and_indel_set);
}

int main(int argc, char
  const * * argv) {

  if (argc != 6) {
    cerr << "USAGE: genesis INPUT.BAM GENE_REGION.BED REFERENCE.FA MASKED_REFERENCE.FA OUTPUT_DIRECTORY\n";
    exit(-1);
  }

  string input_bam = argv[1];
  string gene_region = argv[2];
  string reference_genome = argv[3];
  string masked_reference_genome = argv[4];
  string output_directory = argv[5];

  string input_bam_base_name_with_extension = input_bam.substr(input_bam.find_last_of("\\/") + 1);
  string str(input_bam_base_name_with_extension);

  string input_bam_base_name = str.substr(0, str.find_last_of("."));

  //create output directory
  string create_output_directory = "mkdir -p " + output_directory;
  execute(create_output_directory);

  //extract the base name of the gene_region file
  string gene_region_base_name_with_extension = gene_region.substr(gene_region.find_last_of("\\/") + 1);
  string gene_region_str(gene_region_base_name_with_extension);
  string gene_region_base_name = gene_region_str.substr(0, gene_region_str.find_last_of("."));
/*
  //sort the gene_region
  string sort_gene_region_by_chromosome_and_then_by_position = "sort -k1,1 -k2,2n " + gene_region + " > " + output_directory + "/" + gene_region_base_name + "_sorted.bed";
  execute(sort_gene_region_by_chromosome_and_then_by_position);

  //find and replace all "chrChromosomeNumber" with "ChromosomeNumber" for the sorted gene_region
  string find_and_replace_command = "sed -i 's/chr//g' " + output_directory + "/" + gene_region_base_name + "_sorted.bed";
  execute(find_and_replace_command);
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

  //get the overlap between originally mapped and the gene region
  string get_overlap_between_originally_mapped_and_gene_region = "../ext/bedtools_2.17.0 intersect -abam " + output_directory + "/" + input_bam_base_name + "_mapped-reads.bam -b " + output_directory + "/" + gene_region_base_name + "_sorted.bed > " + output_directory + "/" + input_bam_base_name + "_originally_mapped_overlapped_with_gene_region.bam";
  execute(get_overlap_between_originally_mapped_and_gene_region);

  //generate statistics for the overlap between originally mapped and the gene region
  string generate_statistics_for_overlap_between_originally_mapped_and_gene_region = "../ext/samtools flagstat " + output_directory + "/" + input_bam_base_name + "_originally_mapped_overlapped_with_gene_region.bam > " + output_directory + "/" + input_bam_base_name + "_originally_mapped_overlapped_with_gene_region.stat";
  execute(generate_statistics_for_overlap_between_originally_mapped_and_gene_region);

  //get reads from the originally unmapped bam file
  string get_reads_from_originally_unmapped_bam = "../ext/samtools bam2fq " + output_directory + "/" + input_bam_base_name + "_unmapped-reads.bam > " + output_directory + "/" + input_bam_base_name + "_unmapped-reads.fq";
  execute(get_reads_from_originally_unmapped_bam);
*/
  //map the trimmed reads from the originally unmapped bam file to the masked reference using BWA-MEM
  string map_originally_unampped_read_to_gene_region = "../ext/bwa mem -t 32 " + masked_reference_genome + " " + output_directory + "/" + input_bam_base_name + "_unmapped-reads-trimmed-TruSeq2-SE.fq | ../ext/samtools view -bT " + masked_reference_genome + " - > " + output_directory + "/" + input_bam_base_name + "_bwamem-result-from-unmapped-reads.bam";
  execute(map_originally_unampped_read_to_gene_region);

  //split the newly generated bam file to get the mapped bam file
  string split_newly_mapped_bam_to_mapped_bam = "../ext/samtools view -b -F 0x04 " + output_directory + "/" + input_bam_base_name + "_bwamem-result-from-unmapped-reads.bam > " + output_directory + "/" + input_bam_base_name + "_mapped-from-unmapped.bam";
  execute(split_newly_mapped_bam_to_mapped_bam);

  //generate statistics about the newly mapped reads which were originally unmapped
  string generate_statistics_for_newly_mapped_but_originally_unmapped_reads = "../ext/samtools flagstat " + output_directory + "/" + input_bam_base_name + "_mapped-from-unmapped.bam > " + output_directory + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.stat";
  execute(generate_statistics_for_newly_mapped_but_originally_unmapped_reads);

  //sort the newly mapped bam file
  string sort_newly_mapped_bam = "../ext/samtools sort " + output_directory + "/" + input_bam_base_name + "_mapped-from-unmapped.bam -o " + output_directory + "/" + input_bam_base_name + "_mapped-from-unmapped_sorted.bam";
  execute(sort_newly_mapped_bam);

  //mark PCR duplicates of the newly mapped bam file
  string mark_pcr_duplicate = "java -jar ../ext/MarkDuplicates.jar INPUT=" + output_directory + "/" + input_bam_base_name + "_mapped-from-unmapped_sorted.bam" + " OUTPUT=" + output_directory + "/" + input_bam_base_name + "_mapped-from-unmapped_sorted_dedup.bam" + " METRICS_FILE=" + output_directory + "/" + input_bam_base_name + "_metrics.txt";
  execute(mark_pcr_duplicate);
/*
  //merge the "originally mapped overlapped with the gene region" and "newly mapped" bam files
  string merge_originally_mapped_and_newly_mapped_bam_files = "../ext/samtools merge -f " + output_directory + "/" + input_bam_base_name + "_merged_originally_mapped_overlapped_with_gene_region_and_newly_mapped.bam " + output_directory + "/" + input_bam_base_name + "_originally_mapped_overlapped_with_gene_region.bam " + output_directory + "/" + input_bam_base_name + "_mapped-from-unmapped_sorted_dedup.bam";
  execute(merge_originally_mapped_and_newly_mapped_bam_files);

  //index the merged bam file
  string index_merged_bam_file = "../ext/samtools index " + output_directory + "/" + input_bam_base_name + "_merged_originally_mapped_overlapped_with_gene_region_and_newly_mapped.bam";
  execute(index_merged_bam_file);

  //call snp and indel from the merged bam file using Platypus
  string platypus_command = "python ../ext/Platypus_0.7.9.1/Platypus.py callVariants --bamFiles=" + output_directory + "/" + input_bam_base_name + "_merged_originally_mapped_overlapped_with_gene_region_and_newly_mapped.bam" + " --refFile=" + masked_reference_genome + " --output=" + output_directory + "/" + input_bam_base_name + "_merged_originally_mapped_overlapped_with_gene_region_and_newly_mapped_variants.vcf";
  execute(platypus_command);
*/
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
