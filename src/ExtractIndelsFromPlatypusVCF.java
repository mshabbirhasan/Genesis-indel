
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author shabbir
 */
public class ExtractIndelsFromPlatypusVCF {

    private static void extractIndelFromVCF(String vcfFileName) throws FileNotFoundException, IOException {
        FileReader fileReader = new FileReader(vcfFileName);
        BufferedReader bufferedReader = new BufferedReader(fileReader);

        String outputFileName = vcfFileName.substring(0, vcfFileName.lastIndexOf(".")) + "_indel.vcf";
        FileWriter fileWriter = new FileWriter(outputFileName);

        try {
            String line;

            while ((line = bufferedReader.readLine()) != null) {
                if (line.startsWith("#")) {
                    fileWriter.write(line + "\n");
                } else {
                    String lineArray[] = line.split("\t");
                    int refAlleleLength = lineArray[3].length();
                    if (!lineArray[4].contains(",")) {
                        int altAlleleLength = lineArray[4].length();

                        if (refAlleleLength != altAlleleLength) {
                            fileWriter.write(line + "\n");
                        }
                    } else {
                        String altAlleleArray[] = lineArray[4].split(",");
                        for (int i = 0; i < altAlleleArray.length; i++) {
                            int altAlleLength = altAlleleArray[i].length();
                            if (refAlleleLength != altAlleLength) {
                                fileWriter.write(lineArray[0] + "\t"
                                        + lineArray[1] + "\t"
                                        + lineArray[2] + "\t"
                                        + lineArray[3] + "\t"
                                        + altAlleleArray[i] + "\t"
                                        + lineArray[5] + "\t"
                                        + lineArray[6] + "\t"
                                        + lineArray[7] + "\t"
                                        + lineArray[8] + "\t"
                                        + lineArray[9] + "\n");
                            }
                        }
                    }
                }
            }

            fileReader.close();
            bufferedReader.close();
            fileWriter.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    public static void main(String args[]) {
        if (args.length != 1) {
            System.err.println("USAGE: java ExtractIndelsFromPlatypusVCF VCF_FILE_NAME");
            System.exit(1);
        }

        try {
            extractIndelFromVCF(args[0]);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
