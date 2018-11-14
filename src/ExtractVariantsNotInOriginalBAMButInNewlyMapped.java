
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Shabbir
 */
public class ExtractVariantsNotInOriginalBAMButInNewlyMapped {

  private static void extractNovelVariants(String originallyMappedVcfFileName, String newlyMappedVcfFileName, String outputFileName) throws FileNotFoundException, IOException {
    FileReader fileReader = new FileReader(originallyMappedVcfFileName);
    BufferedReader bufferedReader = new BufferedReader(fileReader);
    
    List<String> originallyMappedIndels = new ArrayList<String>();
    String line;
    
    while((line = bufferedReader.readLine()) != null){
      if(!line.startsWith("#")){
        String lineArray[] = line.split("\t");      
        String chr = lineArray[0];
        String position = lineArray[1];
        String ref = lineArray[3];
        String alt = lineArray[4];
        
        String pattern = chr + "\t" + position + "\t" + ref + "\t" + alt;
        
        originallyMappedIndels.add(pattern);
      }
    }
    
    fileReader.close();
    bufferedReader.close();
    
    fileReader = new FileReader(newlyMappedVcfFileName);
    bufferedReader = new BufferedReader(fileReader);
    
    FileWriter fileWriter = new FileWriter(outputFileName);
    
    while((line = bufferedReader.readLine()) != null){
      if(!line.startsWith("#")){
        String lineArray[] = line.split("\t");      
        String chr = lineArray[0];
        String position = lineArray[1];
        String ref = lineArray[3];
        String alt = lineArray[4];
        
        String pattern = chr + "\t" + position + "\t" + ref + "\t" + alt;
        
        if(!originallyMappedIndels.contains(pattern)){
          fileWriter.write(line + "\n");
        }
      }
    }

    fileReader.close();
    bufferedReader.close();
    fileWriter.close();    
  }

  public static void main(String args[]) {
    if (args.length != 3) {
      System.err.println("USAGE: java ExtractVariantsNotInOriginalBAMButInNewlyMapped VCF_FILE_FOR_ORIGINALLY_MAPPED VCF_FILE_FOR_NEWLY_MAPPED OUTPUT_FILE_NAME");
      System.exit(1);
    }

    try {
      extractNovelVariants(args[0], args[1], args[2]);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
