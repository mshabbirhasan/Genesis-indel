
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
 * @author Shabbir
 */
public class ExtractVariantsWithPassFlagFromAVCF {

  private static void filterVCF(String vcfFileName) throws FileNotFoundException, IOException {
    FileReader fileReader = new FileReader(vcfFileName);
    BufferedReader bufferedReader = new BufferedReader(fileReader);

    FileWriter fileWriter = new FileWriter(vcfFileName.substring(0, vcfFileName.lastIndexOf(".")) + "_PASS_filtered.vcf");
    try {
      String line;
      while ((line = bufferedReader.readLine()) != null) {
        if (line.startsWith("#")) {
          fileWriter.write(line + "\n");
        } else {
          String lineArray[] = line.split("\\t");
          if (lineArray[6].equals("PASS")) {
            fileWriter.write(line + "\n");
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
      System.err.println("USAGE: java ExtractVariantsWithPassFlagFromAVCF VCF_FILE_NAME");
      System.exit(1);
    }
    
    try {
      filterVCF(args[0]);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

}
