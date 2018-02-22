/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import javax.xml.bind.DatatypeConverter;

/**
 *
 * @author worawich
 */
public class FastqUtil {
    
    public static void createIntesectFastqFile(String fq1, String fq2) throws IOException{
        /**
         * This function will find the intersect read sequence between two fastq file
         * The smallest file will be pick as main file store sequence in Map
         * walk through another file and check intersect
         * export intersect read sequence in fastq
         */
        
        File fastq1 = new File(fq1);
        File fastq2 = new File(fq2);
        String exportFile = "";
        String mainFile = "";
        String minorFile = "";
        
        Map<String,ArrayList<String>> mainSeqMap = new LinkedHashMap();
        Map<String,ArrayList<String>> intersectSeqMap = new LinkedHashMap();
        ArrayList<String> mainSeqNameList = new ArrayList();
        
        if(fastq1.length() > fastq2.length()){
            mainFile = fq2;
            minorFile = fq1;
            
            exportFile = fastq2.getParent()+ File.separator + fastq2.getName().split("\\.")[0] + "_Intersect_" + fastq1.getName();
        }else{
            mainFile = fq1;
            minorFile = fq2;
            
            exportFile = fastq1.getParent() + File.separator + fastq1.getName().split("\\.")[0] + "_Intersect_" + fastq2.getName();
        }
        
        Charset charset = Charset.forName("US-ASCII");
        Path mainPath = Paths.get(mainFile);
        Path minorPath = Paths.get(minorFile);
        
        /**
         * Read main fastq File and store in Map<String,String>
         */
        try (BufferedReader reader = Files.newBufferedReader(mainPath, charset)) {
            String line = null;                   
            String name = null;           
            
            int counter = 1;
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    name = null;
                }else{
                    if(line.charAt(0)=='@'&&counter==1){
                        
                        name = line.substring(1);
                        mainSeqNameList.add(name);
                    }else if(counter==2){                    
//                        seq = line;     
                    }else if(counter==3){
//                        strand = line;
                    }else if(counter==4){
//                        quality = line;
//                        
//                        ArrayList<String> fqInfo = new ArrayList();
//                        fqInfo.add(seq);
//                        fqInfo.add(strand);
//                        fqInfo.add(quality);
//                        
//                        mainSeqMap.put(name, fqInfo);
                    }
                }
                
                if(counter==4){
                    counter=1;
                }else{
                    counter++;
                }
            }
        }
        /********************************************************/
        
        /**
         * Read minor fastq File check intersect with mainSeqMap
         * 
         * Store intersect Result
         */
        try (BufferedReader reader = Files.newBufferedReader(minorPath, charset)) {
            String line = null;                   
            String name = null;
            String seq = null;
            String quality = null;
            String strand = null;
            
            boolean intersectFlag = false;
            int counter = 1;
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    name = null;
                }else{
                    if(line.charAt(0)=='@'&&counter==1){
                        
                        name = line.substring(1);
                        
                        if(mainSeqNameList.contains(name)){
                            intersectFlag = true;
                        }else{
                            intersectFlag = false;
                        }

                    }else if(counter==2 && intersectFlag == true){                    
                        seq = line;     
                    }else if(counter==3 && intersectFlag == true){
                        strand = line;
                    }else if(counter==4 && intersectFlag == true){
                        quality = line;
                        
                        ArrayList<String> fqInfo = new ArrayList();
                        fqInfo.add(seq);
                        fqInfo.add(strand);
                        fqInfo.add(quality);
                        
                        intersectSeqMap.put(name, fqInfo);
                    }
                }
                
                if(counter==4){
                    counter=1;
                }else{
                    counter++;
                }    
            }
        }        
        /************************************************************/
        
        /**
         * Export intersect read in fasta file format
         */
        
        FileWriter writer;
        writer = new FileWriter(exportFile);

        for(Map.Entry<String,ArrayList<String>> entry: intersectSeqMap.entrySet()){ 
            ArrayList<String> fqInfo = entry.getValue();
            writer.write("@"+entry.getKey());
            writer.write("\n");
            writer.write(fqInfo.get(0));
            writer.write("\n");
            writer.write(fqInfo.get(1));
            writer.write("\n");
            writer.write(fqInfo.get(2));
            writer.write("\n");
        }
        /**************************************************************/
        
        writer.flush();
        writer.close();
    }
    
    public static void createUnIntesectFastqFile(String fq1, String fq2, String saveFilename) throws IOException, NoSuchAlgorithmException{
        /**
         * This function will find the intersect read sequence between two fastq file
         * fq1 is main file
         * walk through another file and check unIntersect
         * export intersect read sequence in fastq
         * 
         * File can be fq or fq.gz
         * Program will automatically output in fq.gz file if .gz is define in output file name
         */
        
        File fastq1 = new File(fq1);
        File fastq2 = new File(fq2);
        String exportFile = "";
        String mainFile = "";
        String minorFile = "";
        
        Map<String,ArrayList<String>> mainSeqMap = new LinkedHashMap();
        Map<String,ArrayList<String>> unIntersectSeqMap = new LinkedHashMap();
        ArrayList<String> mainSeqNameList = new ArrayList();
        ArrayList<String> mainSeqMd5List = new ArrayList();
        
//        if(fastq1.length() > fastq2.length()){
//            mainFile = fq2;
//            minorFile = fq1;
//            
//            exportFile = fastq2.getParent()+ File.separator + fastq2.getName().split("\\.")[0] + "_unIntersect_" + fastq1.getName();
//        }else{
        mainFile = fq1;     // Be a template for other file to map with (should be the file that did not contain thing that we want)
        minorFile = fq2;    // Map to to template (should contain the thing that we want to extract by un-intersect) the thing in this file that did not apear in template will print out (save as output)

//        exportFile = fastq1.getParent() + File.separator + fastq1.getName().split("\\.")[0] + "_unIntersect_" + fastq2.getName();
//        exportFile = fastq1.getParent() + File.separator + saveFilename;
        exportFile = saveFilename;
//        }
        
        Charset charset = Charset.forName("US-ASCII");
        Path mainPath = Paths.get(mainFile);
        Path minorPath = Paths.get(minorFile);
        
        String[] splitMainFile = mainFile.split("\\.");
        String[] splitMinorFile = minorFile.split("\\.");
        
        
        int numMainSample = 0;
        int numConsiderSample = 0;
        int numOutSample = 0;
        int numCutSample = 0;
        /**
         * Read main fastq File and store in Map<String,String>
         */
        if(splitMainFile[splitMainFile.length-1].equals("gz")){
            FileInputStream fin = new FileInputStream(mainFile);
            GZIPInputStream mainGzip = new GZIPInputStream(fin);
            InputStreamReader mainGzipStream = new InputStreamReader(mainGzip);
            try (BufferedReader reader = new BufferedReader(mainGzipStream)) {
                String line = null;                   
                String name = null;
                String seq = null;
                byte[] seqMd5 = null;

                int counter = 1;
                while ((line = reader.readLine()) != null) {
                    if(line.isEmpty()){
                        name = null;
                    }else{

                        if(line.charAt(0)=='@'&&counter==1){
                            numMainSample++;
                            name = line.substring(1);
                            mainSeqNameList.add(name);
                        }else if(counter==2){                    
                            seq = line;
                            MessageDigest md = MessageDigest.getInstance("MD5");
                            md.update(seq.getBytes());
                            seqMd5 = md.digest();
                            String m = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();
                            if(!mainSeqMd5List.contains(m)){
                                mainSeqMd5List.add(m);
                            }
                        }else if(counter==3){
    //                        strand = line;
                        }else if(counter==4){
    //                        quality = line;
    //                        
    //                        ArrayList<String> fqInfo = new ArrayList();
    //                        fqInfo.add(seq);
    //                        fqInfo.add(strand);
    //                        fqInfo.add(quality);
    //                        
    //                        mainSeqMap.put(name, fqInfo);
                        }
                    }

                    if(counter==4){
                        counter=1;
                    }else{
                        counter++;
                    }
                }
            }
        }else{
            try (BufferedReader reader = Files.newBufferedReader(mainPath, charset)) {
                String line = null;                   
                String name = null;
                String seq = null;
                byte[] seqMd5 = null;

                int counter = 1;
                while ((line = reader.readLine()) != null) {
                    if(line.isEmpty()){
                        name = null;
                    }else{

                        if(line.charAt(0)=='@'&&counter==1){
                            numMainSample++;
                            name = line.substring(1);
                            mainSeqNameList.add(name);
                        }else if(counter==2){                    
                            seq = line;
                            MessageDigest md = MessageDigest.getInstance("MD5");
                            md.update(seq.getBytes());
                            seqMd5 = md.digest();
                            String m = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();
                            if(!mainSeqMd5List.contains(m)){
                                mainSeqMd5List.add(m);
                            }
                        }else if(counter==3){
    //                        strand = line;
                        }else if(counter==4){
    //                        quality = line;
    //                        
    //                        ArrayList<String> fqInfo = new ArrayList();
    //                        fqInfo.add(seq);
    //                        fqInfo.add(strand);
    //                        fqInfo.add(quality);
    //                        
    //                        mainSeqMap.put(name, fqInfo);
                        }
                    }

                    if(counter==4){
                        counter=1;
                    }else{
                        counter++;
                    }
                }
            }
        }
        /********************************************************/
        
        /**
         * Read minor fastq File check intersect with mainSeqMap
         * 
         * Store unintersect Result
         */
        if(splitMinorFile[splitMinorFile.length-1].equals("gz")){
            FileInputStream fin = new FileInputStream(minorFile);
            GZIPInputStream minorGzip = new GZIPInputStream(fin);
            InputStreamReader minorGzipStream = new InputStreamReader(minorGzip);
            try (BufferedReader reader = new BufferedReader(minorGzipStream)) {
                String line = null;                   
                String name = null;
                String seq = null;
                String quality = null;
                String strand = null;
                byte[] seqMd5 = null;

                boolean unIntersectFlag = false;
                int counter = 1;
                while ((line = reader.readLine()) != null) {
                    if(line.isEmpty()){
                        name = null;
                    }else{                    
                        if(line.charAt(0)=='@'&&counter==1){
                            numConsiderSample++;
                            name = line.substring(1);

                            if(!mainSeqNameList.contains(name)){
                                unIntersectFlag = true;
                            }else{
                                unIntersectFlag = false;
                            }

                        }else if(counter==2){                    
                            seq = line;
                            if(unIntersectFlag==true){
                                MessageDigest md = MessageDigest.getInstance("MD5");
                                md.update(seq.getBytes());
                                seqMd5 = md.digest();
                                String Md5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();

                                if(!mainSeqMd5List.contains(Md5)){
                                    unIntersectFlag = true;
                                }else{
                                    unIntersectFlag = false;
                                    numCutSample++;
                                }
                            }
                        }else if(counter==3 && unIntersectFlag == true){
                            strand = line;
                        }else if(counter==4 && unIntersectFlag == true){
                            quality = line;

                            ArrayList<String> fqInfo = new ArrayList();
                            fqInfo.add(seq);
                            fqInfo.add(strand);
                            fqInfo.add(quality);

                            unIntersectSeqMap.put(name, fqInfo);
                            numOutSample++;
                        }
                    }

                    if(counter==4){
                        counter=1;
                    }else{
                        counter++;
                    }    
                }
            }
        }else{
            try (BufferedReader reader = Files.newBufferedReader(minorPath, charset)) {
                String line = null;                   
                String name = null;
                String seq = null;
                String quality = null;
                String strand = null;
                byte[] seqMd5 = null;

                boolean unIntersectFlag = false;
                int counter = 1;
                while ((line = reader.readLine()) != null) {
                    if(line.isEmpty()){
                        name = null;
                    }else{                    
                        if(line.charAt(0)=='@'&&counter==1){
                            numConsiderSample++;
                            name = line.substring(1);

                            if(!mainSeqNameList.contains(name)){
                                unIntersectFlag = true;
                            }else{
                                unIntersectFlag = false;
                            }

                        }else if(counter==2){                    
                            seq = line;
                            if(unIntersectFlag==true){
                                MessageDigest md = MessageDigest.getInstance("MD5");
                                md.update(seq.getBytes());
                                seqMd5 = md.digest();
                                String Md5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();

                                if(!mainSeqMd5List.contains(Md5)){
                                    unIntersectFlag = true;
                                }else{
                                    unIntersectFlag = false;
                                    numCutSample++;
                                }
                            }
                        }else if(counter==3 && unIntersectFlag == true){
                            strand = line;
                        }else if(counter==4 && unIntersectFlag == true){
                            quality = line;

                            ArrayList<String> fqInfo = new ArrayList();
                            fqInfo.add(seq);
                            fqInfo.add(strand);
                            fqInfo.add(quality);

                            unIntersectSeqMap.put(name, fqInfo);
                            numOutSample++;
                        }
                    }

                    if(counter==4){
                        counter=1;
                    }else{
                        counter++;
                    }    
                }
            }
        }        
        /************************************************************/
        
        /**
         * Export intersect read in fasta file format
         */
        String[] splitExportFile = exportFile.split("\\.");
        
        if(splitExportFile[splitExportFile.length-1].equals("gz")){
            FileOutputStream out = new FileOutputStream(exportFile);
            Writer writer = new OutputStreamWriter(new GZIPOutputStream(out),charset);
            for(Map.Entry<String,ArrayList<String>> entry: unIntersectSeqMap.entrySet()){ 
                ArrayList<String> fqInfo = entry.getValue();
                writer.write("@"+entry.getKey());
                writer.write("\n");
                writer.write(fqInfo.get(0));
                writer.write("\n");
                writer.write(fqInfo.get(1));
                writer.write("\n");
                writer.write(fqInfo.get(2));
                writer.write("\n");
            }
            writer.flush();
            writer.close();
        }else{
            FileWriter writer;
            writer = new FileWriter(exportFile);

            for(Map.Entry<String,ArrayList<String>> entry: unIntersectSeqMap.entrySet()){ 
                ArrayList<String> fqInfo = entry.getValue();
                writer.write("@"+entry.getKey());
                writer.write("\n");
                writer.write(fqInfo.get(0));
                writer.write("\n");
                writer.write(fqInfo.get(1));
                writer.write("\n");
                writer.write(fqInfo.get(2));
                writer.write("\n");
            }            
            writer.flush();
            writer.close();
        }
        /**************************************************************/
        System.out.println("num main sample = "+numMainSample);
        System.out.println("num consider sample = "+numConsiderSample);
        System.out.println("num out sample = "+numOutSample);
        System.out.println("num cut sample = "+numCutSample);
    }
    
    
}
