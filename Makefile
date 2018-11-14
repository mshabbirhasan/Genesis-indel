CXX = g++ -std=c++11 

default: all
all: genesis_indel extract_indel_from_vcf extract_novel_variant extract_novel_high_quality_variant

genesis_indel:
	$(CXX) src/genesis-indel.cpp -o bin/genesis-indel -lpthread
	
extract_indel_from_vcf:
	javac -d bin/ src/ExtractIndelsFromPlatypusVCF.java 

extract_novel_variant:
	javac -d bin/ src/ExtractVariantsNotInOriginalBAMButInNewlyMapped.java

extract_novel_high_quality_variant:
	javac -d bin/ src/ExtractVariantsWithPassFlagFromAVCF.java

clean:
	rm -f bin/*

.phony: clean default
