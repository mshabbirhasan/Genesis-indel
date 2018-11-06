#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <thread>
#include <unistd.h>
#include <cstring>

using namespace std;

void execute(string command) {
  int status = system(command.c_str());
  if (status < 0) {
    cerr << "Error: error while executing the command: " << command << endl;
    exit(-1);
  }
}

string createTempDir() {
  char temp_dir_template[] = "./tmpdir.XXXXXX";
  char* temp_dir_name = mkdtemp(temp_dir_template);

  if (temp_dir_name == NULL) {
    cerr << "ERROR: Failed creating the temp directory: " << temp_dir_template << endl;
    exit(-1);
  }

  return temp_dir_name;
}

void remove_temp_dir(string temp_dir_name) {
  if (rmdir(temp_dir_name.c_str()) == -1) {
    cerr << "ERROR: Failed to delete the temp directory: " << temp_dir_name << endl;
    exit(-1);
  }
}

int main(int argc, char* *argv) {
  if (argc != 5) {
    cerr << "USAGE: genesis-indel INPUT.BAM REFERENCE.FA OUTPUT_FILENAME -reportSNP=true\n";
    exit(-1);
  }

  string input_bam = argv[1];
  string reference_genome = argv[2];
  string output_file_name = argv[3];
  string reportSNP = argv[4];

  string input_bam_base_name_with_extension = input_bam.substr(input_bam.find_last_of("\\/") + 1);
  string str(input_bam_base_name_with_extension);

  string input_bam_base_name = str.substr(0, str.find_last_of("."));

  // create a temp directory
  string temp_dir_name = createTempDir();

  // split the input BAM file into mapped and unmapped bam file
  string split_to_mapped_bam_command = "../ext/samtools view -b -F 0x04 " + input_bam + " > "
          + temp_dir_name + "/" + input_bam_base_name + "_mapped-reads.bam";
  string split_to_unmapped_bam_command = "../ext/samtools view -b -f 0x04 " + input_bam + " > "
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.bam";

  thread split_to_mapped_bam(execute, split_to_mapped_bam_command);
  thread split_to_unmapped_bam(execute, split_to_unmapped_bam_command);

  split_to_mapped_bam.join(); // pauses until split_to_mapped_bam finishes
  split_to_unmapped_bam.join(); // pauses until split_to_unmapped_bam finishes  

  // get reads from the originally unmapped bam file
  string get_reads_from_originally_unmapped_bam = "../ext/samtools bam2fq " + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.bam > "
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.fq";
  execute(get_reads_from_originally_unmapped_bam);
  
  // do the quality control of the unmapped read using Trimmomatic
  string trim_unmapped_read_using_trimmomatic_command = "java -jar ../ext/trimmomatic.jar SE -phred33 " 
          + temp_dir_name + "/" + input_bam_base_name + "_unmapped-reads.fq " 
          + temp_dir_name + "/" + input_bam_base_name 
          + "_unmapped-reads-trimmed-TruSeq2-SE.fq ILLUMINACLIP:../ext/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
  execute(trim_unmapped_read_using_trimmomatic_command);
  
  // remove temp directory
  remove_temp_dir(temp_dir_name);
}
