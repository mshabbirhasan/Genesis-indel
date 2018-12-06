# Genesis-indel
Genesis-indel is a computational pipeline to explore the unmapped reads to identify mutations that are initially missed in the original alignment. The genes containing such mutations can be investigated to gain important biological insights.

**Author: Mohammad Shabbir Hasan**

**Ph.D. Student, Department of Computer Science**

**Virginia Tech, Blacksburg, VA 24060, USA.**

**Email: shabbir5@vt.edu**

## System Requirement
- gcc: version 4.9 or above
- jdk: version 1.8 or above
- python: Python 2.6 or above
- bedtools: version 2.27 or above (Can be downloaded from https://bedtools.readthedocs.io/en/latest/content/installation.html)
- blat: version 36 X 2 or above (Executables can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/)
- bwa: version 0.7.15 or above (Can be downloaded from https://github.com/lh3/bwa)
- samtools: version 1.4 or above (Can be downloaded from http://www.htslib.org/)

**Precompiled executables are provided in the ext/ directory that are tested to work with Ubuntu 14.04 LTS.**

## Build and install
### Preparation:
- After you download and install bedtools, add bedtools to $PATH variable using the command 
    ```PATH=path_to_bedtools_bin:$PATH```.
    **Alternatively**, if you use the precompiled bedtools executable provided in ext/, you can add it to the $PATH variable as follows:
        ```PATH=ext/:$PATH```. This will add other precompiled executables provided in ext/ to the $PATH variable.
- Add blat, bwa, and samtools to the $PATH variable in the same way.
### Build and install        
**Step 1:** 
    ```make clean``` // This will remove previous binaries.
    
**Step 2:**
    ``` make ``` // This will create binaries in the bin/ directory.
    
**Step 3:**
    ``` make install ``` // This will give executable permission to the external binaries.

## Run Genesis-indel
### Input: 
**1. Alignment file in BAM format.**
> There should be a corresponding index file (.bai file) of the input BAM in the same direcory. You can create an index file using the following command:

``` samtools index input.bam ``` 

***Note:*** samtools is included in the ext/ directory.

**2. Reference genome in FASTA format.**
> There should be corresponding index files (.fai, .amb, .ann, .bwt, .pac, .sa) of the input FASTA in the same directory. You can create the index files using the follwoing commands:

``` 
samtools faidx reference.fasta
bwa index reference.fasta 
```

***Note:*** You need to create the index for reference only once.

***Note:*** bwa is included in the ext/ directory.

### Output: 
**A VCF file containing the novel high-quality variants.**

## Example
``` ./genesis-indel input.bam reference.fasta output.vcf -reportSNP=false ```

This command will report novel high-quality indels only.

``` ./genesis-indel input.bam reference.fasta output.vcf -reportSNP=true ```

This command will report both novel high-quality SNPs and indels.

## Citing Genesis-indel
If you use Genesis-indel, please cite the following paper:

Mohammad Shabbir Hasan, Xiaowei Wu, and Liqing Zhang. "Uncovering missed indels by leveraging unmapped reads." TBD.