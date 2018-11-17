#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <thread>
#include <unistd.h>
#include <cstring>
#include <algorithm>
#include <dirent.h>
#include <vector>
#include <memory>

using namespace std;

void execute(string command) {
  int status = system(command.c_str());
  if (status < 0) {
    cerr << "Error: error while executing the command " << command << strerror(errno) << endl;
    exit(-1);
  }
}

string createTempDir() {
  char temp_dir_template[] = "./tmpdir.XXXXXX";
  char* temp_dir_name = mkdtemp(temp_dir_template);

  if (temp_dir_name == NULL) {
    cerr << "ERROR: Failed creating the temp directory " << temp_dir_template << " : " << strerror(errno) << endl;
    exit(-1);
  }
  return temp_dir_name;
}

vector<string> get_file_list(string directory_name) {
  vector<string> files;
  shared_ptr<DIR> directory_ptr(opendir(directory_name.c_str()), [](DIR * directory_name) {
    directory_name && closedir(directory_name); });
  struct dirent *dirent_ptr;
  if (!directory_ptr) {
    cout << "ERROR: Can not open directory " << directory_name << " : " << strerror(errno) << endl;
    exit(-1); // returns an empty vector
  }

  while ((dirent_ptr = readdir(directory_ptr.get())) != nullptr) {
    string file_name = string(dirent_ptr -> d_name);
    if (file_name != "." && file_name != "..") {
      files.push_back(file_name);
    }
  }
  return files;
}

void remove_temp_dir(string temp_dir_name) {
  vector<string> files = get_file_list(temp_dir_name);
  for (const auto& file : files) {
    string file_to_delete = temp_dir_name + "/" + file;
    if (remove(file_to_delete.c_str()) != 0) {
      cerr << "ERROR: Failed to delete the temp file " << file_to_delete << " : " << strerror(errno) << endl;
      exit(-1);
    }
  }

  if (rmdir(temp_dir_name.c_str()) == -1) {
    cerr << "ERROR: Failed to delete the temp directory " << temp_dir_name << " : " << strerror(errno) << endl;
    exit(-1);
  }
}

void copy_output(string srcFile, string dstFile) {
  ifstream src(srcFile, ios::binary);
  ofstream dst(dstFile, ios::binary);

  if (!(dst << src.rdbuf())) {
    cerr << "ERROR: Failed to copy output " << srcFile << " to the output file " << dstFile << " : " << strerror(errno) << endl;
    src.close();
    dst.close();
    exit(-1);
  }
  src.close();
  dst.close();
}

int main(int argc, char* *argv) {
  if (argc != 5) {
    cerr << "USAGE: genesis-indel INPUT.BAM REFERENCE.FA OUTPUT_FILENAME -reportSNP=false\n";
    exit(-1);
  }

  string input_bam = argv[1];
  string reference_genome = argv[2];
  string output_file_name = argv[3];
  string reportSNP = argv[4];
  string flag = reportSNP.substr(11, reportSNP.length());
  transform(flag.begin(), flag.end(), flag.begin(), ::tolower);

  string input_bam_base_name_with_extension = input_bam.substr(input_bam.find_last_of("\\/") + 1);
  string str(input_bam_base_name_with_extension);

  string input_bam_base_name = str.substr(0, str.find_last_of("."));

  // create a temp directory
  string temp_dir_name = createTempDir();

  // extract the unmapped bam from the input bam
  string command_to_extract_unmapped_bam_from_input_bam = "../ext/samtools view -b -f 0x04 " + input_bam + " > "
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.bam";
  execute(command_to_extract_unmapped_bam_from_input_bam);

  // get reads from the originally unmapped bam file
  string command_to_extract_reads_from_originally_unmapped_bam = "../ext/samtools bam2fq "
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.bam > "
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.fq";
  execute(command_to_extract_reads_from_originally_unmapped_bam);

  // do the quality control of the unmapped read using Trimmomatic
  string command_to_trim_unmapped_read_using_trimmomatic = "java -jar ../ext/trimmomatic.jar SE -phred33 "
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.fq "
          + temp_dir_name + "/" + input_bam_base_name
          + "_unmapped-reads-trimmed-TruSeq2-SE.fq ILLUMINACLIP:../ext/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
  execute(command_to_trim_unmapped_read_using_trimmomatic);

  // map the trimmed reads from the originally unmapped reads to the reference using BWA-MEM
  string command_to_map_originally_unampped_read_to_reference_genome = "../ext/bwa mem -t 32 " + reference_genome + " "
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads-trimmed-TruSeq2-SE.fq | ../ext/samtools view -bT "
          + reference_genome + " - > " + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.bam";
  execute(command_to_map_originally_unampped_read_to_reference_genome);

  // extract reads that are still unmapped after BWA-MEM
  string command_to_extract_still_unmapped_reads_after_BWA_MEM = "../ext/samtools view -b -f 0x04 "
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.bam > "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM.bam";
  execute(command_to_extract_still_unmapped_reads_after_BWA_MEM);

  // convert still unmapped reads to fasta
  string command_to_convert_still_unammped_reads_to_fasta = "../ext/samtools fasta "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM.bam > "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM.fasta";
  execute(command_to_convert_still_unammped_reads_to_fasta);

  // run BLAT for still unmapped reads
  string command_to_align_still_unmapped_reads_using_BLAT = "../ext/blat " + reference_genome + " "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM.fasta "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM_aligned_by_BLAT.psl";
  execute(command_to_align_still_unmapped_reads_using_BLAT);

  // convert PSL to BED
  string command_to_convert_psl_to_bed = "../ext/pslToBed "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM_aligned_by_BLAT.psl "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM_aligned_by_BLAT.bed";
  execute(command_to_convert_psl_to_bed);

  // convert BED to BAM
  string command_to_get_the_chromosome_size = "awk -v OFS='\t' {'print $1, $2'} " + reference_genome + ".fai > " + temp_dir_name + "/chr_size.txt";
  execute(command_to_get_the_chromosome_size);
  string command_to_convert_bed_to_bam = "../ext/bedtools bedtobam -bed12 -i "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM_aligned_by_BLAT.bed -g "
          + temp_dir_name + "/chr_size.txt > " + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM_aligned_by_BLAT.bam";
  execute(command_to_convert_bed_to_bam);

  // merge BAM from BWA-MEM and BAM from BLAT
  string command_to_merge_bwa_mem_bam_and_blat_mem = "../ext/samtools merge -f "
          + temp_dir_name + "/" + input_bam_base_name + "_BWA-MEM_and_BLAT_merged.bam "
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.bam "
          + temp_dir_name + "/" + input_bam_base_name + "_still_unmapped_after_BWA-MEM_aligned_by_BLAT.bam";
  execute(command_to_merge_bwa_mem_bam_and_blat_mem);

  // sort the newly mapped bam file
  string command_to_sort_newly_mapped_bam = "../ext/samtools sort -@ 4 "
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped.bam -o "
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted.bam";
  execute(command_to_sort_newly_mapped_bam);

  // mark PCR duplicates of the newly mapped bam file
  string command_to_mark_pcr_duplicate = "java -jar ../ext/MarkDuplicates.jar INPUT="
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted.bam" + " OUTPUT="
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup.bam" + " METRICS_FILE="
          + temp_dir_name + "/" + input_bam_base_name + "_metrics.txt";
  execute(command_to_mark_pcr_duplicate);

  // index the originally unmapped newly mapped sorted dedup bam file
  string command_to_index_merged_bam_file = "../ext/samtools index "
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup.bam";
  execute(command_to_index_merged_bam_file);

  // call SNPs and Indel from the originally mapped reads and also from originally unmapped newly mapped sorted dedup bam file using Platypus
  string command_to_call_variant_from_original_bam = "python ../ext/platypus/Platypus.py callVariants --bamFiles="
          + input_bam + " --refFile=" + reference_genome + " --output="
          + temp_dir_name + "/" + input_bam_base_name + "_original_bam_variant.vcf";

  string command_to_call_variant_from_new_bam = "python ../ext/platypus/Platypus.py callVariants --bamFiles="
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup.bam" + " --refFile=" + reference_genome + " --output="
          + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup_variants.vcf";

  thread call_variant_from_original_bam(execute, command_to_call_variant_from_original_bam);
  thread call_variant_from_new_bam(execute, command_to_call_variant_from_new_bam);

  call_variant_from_original_bam.join();
  call_variant_from_new_bam.join();

  if (flag == "false") {
    // extract indel from Platypus VCF for originally mapped reads and also from newly mapped reads
    string command_to_extract_indel_from_variants_of_original_bam = "java ExtractIndelsFromPlatypusVCF "
            + temp_dir_name + "/" + input_bam_base_name + "_original_bam_variant.vcf";

    string command_to_extract_indel_from_variants_of_new_bam = "java ExtractIndelsFromPlatypusVCF "
            + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup_variants.vcf";

    thread extract_variant_from_original_bam(execute, command_to_extract_indel_from_variants_of_original_bam);
    thread extract_variant_from_new_bam(execute, command_to_extract_indel_from_variants_of_new_bam);

    extract_variant_from_original_bam.join();
    extract_variant_from_new_bam.join();

    // extract novel indels
    string command_to_extract_novel_indel = "java ExtractVariantsNotInOriginalBAMButInNewlyMapped "
            + temp_dir_name + "/" + input_bam_base_name + "_original_bam_variant_indel.vcf "
            + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup_variants_indel.vcf "
            + temp_dir_name + "/" + input_bam_base_name + "_novel_indel.vcf";
    execute(command_to_extract_novel_indel);

    // extract novel high-quality indels are report it in the output
    string command_to_extract_novel_high_quality_indel = "java ExtractVariantsWithPassFlagFromAVCF "
            + temp_dir_name + "/" + input_bam_base_name + "_novel_indel.vcf";
    execute(command_to_extract_novel_high_quality_indel);

    // copy the novel high-quality indels to output
    string novel_high_quality_indel_file_name = temp_dir_name + "/" + input_bam_base_name + "_novel_indel_PASS_filtered.vcf";
    copy_output(novel_high_quality_indel_file_name, output_file_name);
  } else if (flag == "true") {
    // extract novel SNPs and indels
    string command_to_extract_novel_snp_and_indel = "java ExtractVariantsNotInOriginalBAMButInNewlyMapped "
            + temp_dir_name + "/" + input_bam_base_name + "_original_bam_variant.vcf "
            + temp_dir_name + "/" + input_bam_base_name + "_originally_unmapped_newly_mapped_sorted_dedup_variants.vcf "
            + temp_dir_name + "/" + input_bam_base_name + "_novel_SNP_and_indel.vcf";
    execute(command_to_extract_novel_snp_and_indel);

    // extract novel high-quality SNPs and indels
    string command_to_extract_novel_high_quality_snp_and_indel = "java ExtractVariantsWithPassFlagFromAVCF "
            + temp_dir_name + "/" + input_bam_base_name + "_novel_SNP_and_indel.vcf";
    execute(command_to_extract_novel_high_quality_snp_and_indel);

    // copy the novel high-quality indels to output
    string novel_high_quality_SNP_and_indel_file_name = temp_dir_name + "/" + input_bam_base_name + "_novel_SNP_and_indel_PASS_filtered.vcf";
    copy_output(novel_high_quality_SNP_and_indel_file_name, output_file_name);
  }

  // remove temp directory
  remove_temp_dir(temp_dir_name);
}
