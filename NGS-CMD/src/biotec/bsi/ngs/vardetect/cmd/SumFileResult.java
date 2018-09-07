/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.VariationV2;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */
public class SumFileResult {
    public static void main (String args[]) throws IOException{
        /**
         * create sum file of all sv type
         */
        Charset charset = Charset.forName("US-ASCII");
        Boolean firstFlag = true;
        
        /************/
        
        /**
         * Modified map to new file name
         */
        Map<String,String> nameMap = new HashMap();
        String name_file = "/Volumes/4TB_WD/TB/Martin_new_list/TB_remove_dup_bam/Sample_fullName.txt";
        File nameFile = new File(name_file);
        Path pathF = Paths.get(nameFile.getPath());
        try (BufferedReader reader = Files.newBufferedReader(pathF, charset)) {
            String line = null;    
            while ((line = reader.readLine()) != null) {
                String[] portion = line.split("_");
                nameMap.put(portion[0], line);
            }
        }
        
        /**************/
        
        /**
         * start loop through dir
         */
        String myDirectoryPath = "/Volumes/4TB_WD/TB/Martin_new_list/TB_remove_dup_bam/TB_res";
        
        FileWriter writer;

        File dir = new File(myDirectoryPath);
        File[] directoryListing = dir.listFiles();
        if (directoryListing != null) {
            // loop in initial path dir
            for (File child : directoryListing) {
                // initial path folder level
                String sampleName = child.getName();
                /**
                * Modified map to new file name
                */
                sampleName = nameMap.get(sampleName);
                /**************/               
                File nextDir = new File(child.getPath());
                File[] nextDirList = nextDir.listFiles();
                if (nextDirList != null) {
                    // loop in sample dir
                    for (File nextChild : nextDirList) {
                        // sample folder level
                        String svType = nextChild.getName();
                        File svDir = new File(nextChild.getPath());
                        File[] svDirList = svDir.listFiles();
                        if(svDirList != null){
                            // loop in sv dir
                            for (File csvFile : svDirList) {
                                // sv folder level
                                if(csvFile.getName().split("\\.")[1].equals("csv")){
                                    //read csv file
                                    Path path = Paths.get(csvFile.getPath());
                                    
                                    File sumFile = new File(myDirectoryPath + File.separator + "sumResult_"+svType+".csv"); //File object        
                                    if(sumFile.exists()){
//                                        ps = new PrintStream(new FileOutputStream(filename,true));
                                        writer = new FileWriter(sumFile,true);
                                        firstFlag = false;
//                                        writer = new FileWriter(summaryReport); //not append
                                    }else{
//                                        ps = new PrintStream(filename);
                                        writer = new FileWriter(sumFile);
                                        firstFlag = true;
                                    }
                                    
                                    
                                    
                                    try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
                                        String line = null;    
                                        int count = 0;

//                                        System.out.println("reading csv File : "+csvFile.getName());
                                        while ((line = reader.readLine()) != null) {
                                            
                                            if(firstFlag == true && count == 0){
                                                
                                                String modifLine = "SampleName,"+line;
                                                writer.write(modifLine);
                                                writer.write("\n");
                                                
                                            }else if(count > 0){
                                                String modifLine = sampleName+","+line;
                                                writer.write(modifLine);
                                                writer.write("\n");
                                            }
                                            
                                            count++;  
                                        }
                                    }
                                    writer.flush();
                                    writer.close();
                                }
                            }
                        }
                    }
                }
                else{

                }
             }
        } else {
          // Handle the case where dir is not really a directory.
          // Checking dir.isDirectory() above would not be sufficient
          // to avoid race conditions with another process that deletes
          // directories.
        } 
    }
}
