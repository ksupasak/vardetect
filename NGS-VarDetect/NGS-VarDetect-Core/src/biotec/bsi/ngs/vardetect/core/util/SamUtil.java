/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */
public class SamUtil {
    
    public static void createIntesectSamFile(String sam1, String sam2) throws IOException{
        /**
         * This function will find the intersect read sequence between two sam file
         * The smallest file will be pick as main file store sequence in Map
         * walk through another file and chack intersect
         * export intersect read sequence in sam
         */
        
        File samFile1 = new File(sam1);
        File samFile2 = new File(sam2);
        String exportFile = "";
        String mainFile = "";
        String minorFile = "";
        
        ArrayList<String> mainSeq = new ArrayList();
        ArrayList<String> intersectSeq = new ArrayList();
        ArrayList<String> header = new ArrayList();
  
        if(samFile1.length() > samFile2.length()){
            mainFile = sam2;
            minorFile = sam1;
            
            exportFile = samFile2.getParent()+ File.separator + samFile2.getName().split("\\.")[0] + "_Intersect_" + samFile1.getName();
        }else{
            mainFile = sam1;
            minorFile = sam2;
            
            exportFile = samFile1.getParent() + File.separator + samFile1.getName().split("\\.")[0] + "_Intersect_" + samFile2.getName();
        }
        
        Charset charset = Charset.forName("US-ASCII");
        Path mainPath = Paths.get(mainFile);
        Path minorPath = Paths.get(minorFile);
        
        /**
         * Read main fasta File and store in Map<String,String>
         */
        try (BufferedReader reader = Files.newBufferedReader(mainPath, charset)) {
            String line = null;                   
            String name = null;
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)!='@'){                        
                        name = line.split("\t")[0];
                        mainSeq.add(name);
                    }else{
                        header.add(line);
                    }
                }    
            }
        }
        /********************************************************/
        
        /**
         * Read minor fasta File check intersect with mainSeqMap
         * 
         * Store intersect Result
         */
        try (BufferedReader reader = Files.newBufferedReader(minorPath, charset)) {
            String line = null;                   
            String name = null;            
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)!='@'){                        
                        name = line.split("\t")[0];
                        if(mainSeq.contains(name)){
                            intersectSeq.add(line);
                        }
                    }
                }    
            }
        }        
        /************************************************************/
        
        /**
         * Export intersect read in sam file format
         */
        
        FileWriter writer;
        writer = new FileWriter(exportFile);
        
        /**
         * write header
         */
        for(String headerData : header){
            writer.write(headerData);
            writer.write("\n");
        }
        
        /**
         * write data
         */
        for(String data : intersectSeq){
            writer.write(data);
            writer.write("\n");
        }
        
        writer.flush();
        writer.close();
        /**************************************************************/
    }
    
    public static void createUnIntesectSamFile(String sam1, String sam2) throws IOException{
        /**
         * This function will find the intersect read sequence between two fasta file
         * The smallest file will be pick as main file store sequence in Map
         * walk through another file and chack intersect
         * export intersect read sequence in fasta
         */
        
        File samFile1 = new File(sam1);
        File samFile2 = new File(sam2);
        String exportFile = "";
        String mainFile = "";
        String minorFile = "";
        
        ArrayList<String> mainSeq = new ArrayList();
        ArrayList<String> unIntersectSeq = new ArrayList();
        ArrayList<String> header = new ArrayList();
  
        if(samFile1.length() > samFile2.length()){
            mainFile = sam2;
            minorFile = sam1;
            
            exportFile = samFile2.getParent()+ File.separator + samFile2.getName().split("\\.")[0] + "_unIntersect_" + samFile1.getName();
        }else{
            mainFile = sam1;
            minorFile = sam2;
            
            exportFile = samFile1.getParent() + File.separator + samFile1.getName().split("\\.")[0] + "_unIntersect_" + samFile2.getName();
        }
        
        Charset charset = Charset.forName("US-ASCII");
        Path mainPath = Paths.get(mainFile);
        Path minorPath = Paths.get(minorFile);
        
        /**
         * Read main fasta File and store in Map<String,String>
         */
        try (BufferedReader reader = Files.newBufferedReader(mainPath, charset)) {
            String line = null;                   
            String name = null;
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)!='@'){                        
                        name = line.split("\t")[0];
                        mainSeq.add(name);
                    }else{
                        header.add(line);
                    }
                }    
            }
        }
        /********************************************************/
        
        /**
         * Read minor fasta File check intersect with mainSeqMap
         * 
         * Store intersect Result
         */
        try (BufferedReader reader = Files.newBufferedReader(minorPath, charset)) {
            String line = null;                   
            String name = null;            
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)!='@'){                        
                        name = line.split("\t")[0];
                        if(!mainSeq.contains(name)){
                            unIntersectSeq.add(line);
                        }
                    }
                }    
            }
        }        
        /************************************************************/
        
        /**
         * Export intersect read in sam file format
         */
        
        FileWriter writer;
        writer = new FileWriter(exportFile);
        
        /**
         * write header
         */
        for(String headerData : header){
            writer.write(headerData);
            writer.write("\n");
        }
        
        /**
         * write data
         */
        for(String data : unIntersectSeq){
            writer.write(data);
            writer.write("\n");
        }
        
        writer.flush();
        writer.close();
        /**************************************************************/
    }
}
