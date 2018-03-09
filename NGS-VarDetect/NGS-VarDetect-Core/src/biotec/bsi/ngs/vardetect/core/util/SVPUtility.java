/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.CombineReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ThreadPool;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import javax.xml.bind.DatatypeConverter;

/**
 *
 * @author worawich
 */
public class SVPUtility {
    
    public static void modifyFile(String file) throws IOException{
        /**
         * read file and modify then write to new file
         */
        
        FileWriter writer;
        String[] dummysaveFile = file.split("\\.");
        String saveFile = dummysaveFile[0] + "_Modif." + dummysaveFile[1];
        /**
         * Check File existing
         */
        
        File f = new File(saveFile); //File object        
        if(f.exists()){
            // append if exist
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(saveFile,true);
        }else{
            // create new
//            ps = new PrintStream(filename);
            writer = new FileWriter(saveFile);
        }
        writer = new FileWriter(saveFile);
        
        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(file);
        
        StringBuffer seq = new StringBuffer();

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;
            String[] data = null; 
            while ((line = reader.readLine()) != null) {
                data = line.split("_");
                writer.write(data[0]);
                writer.write("\n");
            }
        }
        
        writer.flush();
        writer.close();
    }
    
    public static void addUnmapFlagToBam(String inputFile,String targetFile,String fileFormat) throws NoSuchAlgorithmException, IOException{
        /**
         * Force add unMap bit flag and * at cigar to bam file
         * user specify file format in sam or bam
         */
        int count = 0;
        String outputFile = inputFile.split("\\.")[0] + "_modif."+fileFormat;
        File bamFile = new File(inputFile);
        String[] targetSTRComp = targetFile.split("\\.");
        String targetFormat = targetSTRComp[targetSTRComp.length-1];
        
        ArrayList<String> targetSeqMd5List = new ArrayList();
        Map<String,Boolean> targetSeqMd5Map = new LinkedHashMap();
            
        if(targetFormat.equals("bam")){
            File targetBamFile = new File(targetFile);
            SamReaderFactory samReaderFactoryInput = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
            SamReader samReaderInput = samReaderFactoryInput.open(bamFile);
            SamReaderFactory samReaderFactoryTarget = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
            SamReader samReaderTarget = samReaderFactoryTarget.open(targetBamFile);

            File outputBamFile = new File(outputFile);
            SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReaderInput.getFileHeader(),true, outputBamFile);

            /**
             * loop read target bam file
             * create sequence md5 for matching purpose
             */
            
            byte[] seqMd5 = null;          
            int counter = 1;
            SAMRecordIterator itTarget = samReaderTarget.iterator();

            while (itTarget.hasNext()){
                SAMRecord data = itTarget.next();
                String seq = data.getReadString();
                MessageDigest md = MessageDigest.getInstance("MD5");
                md.update(seq.getBytes());
                seqMd5 = md.digest();
                String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();                

                if(!targetSeqMd5Map.containsKey(strSeqMd5)){
                    targetSeqMd5Map.put(strSeqMd5,true);
                }         
            }
            /**********************************************/

            /**
             * loop main bam file 
             * Force add 1 at bit 3 of sam flags (only on target seq)
             */
            SAMRecordIterator itMain = samReaderInput.iterator();

            while(itMain.hasNext()){
                SAMRecord data = itMain.next();

                String seq = data.getReadString();
                MessageDigest md = MessageDigest.getInstance("MD5");
                md.update(seq.getBytes());
                seqMd5 = md.digest();
                String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();                

                if(targetSeqMd5Map.containsKey(strSeqMd5)){
                    int numFlag = data.getFlags();
                    if(numFlag == 3 || numFlag==2){
                        StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.append('1');
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }else if(numFlag==1 || numFlag==0){
                        StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.append("10");
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }else{
                        StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.setCharAt((strFlag.length()-1)-2,'1');
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }
                    count++;
                }
                outputSam.addAlignment(data);
            }
            /*****************************************************/

            outputSam.close();
            samReaderInput.close();
            samReaderTarget.close();

            System.out.println("num read modified = "+count);
        }else if(targetFormat.equals("gz")&&targetSTRComp[targetSTRComp.length-2].equals("fq")){
            
            SamReaderFactory samReaderFactoryInput = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
            SamReader samReaderInput = samReaderFactoryInput.open(bamFile);
            
            FileInputStream fin = new FileInputStream(targetFile);
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
                        if(counter==2){                    
                            seq = line;
                            MessageDigest md = MessageDigest.getInstance("MD5");
                            md.update(seq.getBytes());
                            seqMd5 = md.digest();
                            String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();
                            if(!targetSeqMd5Map.containsKey(strSeqMd5)){
                                targetSeqMd5Map.put(strSeqMd5,true);
                            }
                        }
                    }
                    if(counter==4){
                        counter=1;
                    }else{
                        counter++;
                    }
                }
            }
            

            File outputBamFile = new File(outputFile);
            SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReaderInput.getFileHeader(),true, outputBamFile);
            /**
             * loop main bam file 
             * Force add 1 at bit 3 of sam flags (only on target seq)
             */
            SAMRecordIterator itMain = samReaderInput.iterator();

            while(itMain.hasNext()){
                byte[] seqMd5 = null;
                SAMRecord data = itMain.next();
                StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(data.getFlags()));
                boolean reverseFlag = false;
                if(strFlag.charAt((strFlag.length()-1)-4)=='1'){
                    reverseFlag = true;
                }
                
                String seq = data.getReadString();
                
                if(reverseFlag==true){
                    String invSeq = SequenceUtil.inverseSequence(seq);                                // Do invert sequence (ATCG => GCTA)
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)
                    seq = compSeq;
                }
                
                MessageDigest md = MessageDigest.getInstance("MD5");
                md.update(seq.getBytes());
                seqMd5 = md.digest();
                String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();                

                if(targetSeqMd5Map.containsKey(strSeqMd5)){
                    int numFlag = data.getFlags();
                    if(numFlag == 3 || numFlag==2){
//                        strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.append('1');
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }else if(numFlag==1 || numFlag==0){
//                        strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.append("10");
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }else{
//                        strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.setCharAt((strFlag.length()-1)-2,'1');
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }
                    count++;
                }
                outputSam.addAlignment(data);
            }
            /*****************************************************/

            outputSam.close();
            samReaderInput.close();

            System.out.println("num read modified = "+count);
            
        }else if(targetFormat.equals("fq")){
            SamReaderFactory samReaderFactoryInput = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
            SamReader samReaderInput = samReaderFactoryInput.open(bamFile);
            
            Path mainPath = Paths.get(targetFile);
            Charset charset = Charset.forName("US-ASCII");
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
                        if(counter==2){                    
                            seq = line;
                            MessageDigest md = MessageDigest.getInstance("MD5");
                            md.update(seq.getBytes());
                            seqMd5 = md.digest();
                            String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();
                            if(!targetSeqMd5Map.containsKey(strSeqMd5)){
                                targetSeqMd5Map.put(strSeqMd5,true);
                            }
                        }
                    }

                    if(counter==4){
                        counter=1;
                    }else{
                        counter++;
                    }
                }
            }
            
            File outputBamFile = new File(outputFile);
            SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReaderInput.getFileHeader(),true, outputBamFile);
            /**
             * loop main bam file 
             * Force add 1 at bit 3 of sam flags (only on target seq)
             */
            SAMRecordIterator itMain = samReaderInput.iterator();

            while(itMain.hasNext()){
                byte[] seqMd5 = null;
                SAMRecord data = itMain.next();
                StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(data.getFlags()));
                boolean reverseFlag = false;
                if(strFlag.charAt((strFlag.length()-1)-4)=='1'){
                    reverseFlag = true;
                }
                
                String seq = data.getReadString();
                
                if(reverseFlag==true){
                    String invSeq = SequenceUtil.inverseSequence(seq);                                // Do invert sequence (ATCG => GCTA)
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)
                    seq = compSeq;
                }
                
                MessageDigest md = MessageDigest.getInstance("MD5");
                md.update(seq.getBytes());
                seqMd5 = md.digest();
                String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();                

                if(targetSeqMd5Map.containsKey(strSeqMd5)){
                    int numFlag = data.getFlags();
                    if(numFlag == 3 || numFlag==2){
//                        strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.append('1');
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }else if(numFlag==1 || numFlag==0){
//                        strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.append("10");
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }else{
//                        strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                        strFlag.setCharAt((strFlag.length()-1)-2,'1');
                        data.setFlags(Integer.parseInt(strFlag.toString(),2));
                    }
                    count++;
                }
                outputSam.addAlignment(data);
            }
            /*****************************************************/

            outputSam.close();
            samReaderInput.close();

            System.out.println("num read modified = "+count);
            
        }
    }
    
    public static void addUnmapFlagToBamOnlyInputBam(String inputFile,String targetFile,String fileFormat) throws NoSuchAlgorithmException, IOException{
        /**
         * Force add unMap bit flag and * at cigar to bam file
         * user specify file format in sam or bam
         */
        int count = 0;
        String outputFile = inputFile.split("\\.")[0] + "_modif."+fileFormat;
        File bamFile = new File(inputFile);
        File targetBamFile = new File(targetFile);
        SamReaderFactory samReaderFactoryInput = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader samReaderInput = samReaderFactoryInput.open(bamFile);
        SamReaderFactory samReaderFactoryTarget = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader samReaderTarget = samReaderFactoryTarget.open(targetBamFile);
        
        File outputBamFile = new File(outputFile);
        SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReaderInput.getFileHeader(),true, outputBamFile);
        
        /**
         * loop read target bam file
         * create sequence md5 for matching purpose
         */
        ArrayList<String> targetSeqMd5List = new ArrayList();
        Map<String,Boolean> targetSeqMd5Map = new LinkedHashMap();
        byte[] seqMd5 = null;          
        int counter = 1;
        SAMRecordIterator itTarget = samReaderTarget.iterator();
        
        while (itTarget.hasNext()){
            SAMRecord data = itTarget.next();
            String seq = data.getReadString();
            MessageDigest md = MessageDigest.getInstance("MD5");
            md.update(seq.getBytes());
            seqMd5 = md.digest();
            String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();                
            
            if(!targetSeqMd5Map.containsKey(strSeqMd5)){
                targetSeqMd5Map.put(strSeqMd5,true);
            }         
        }
        /**********************************************/
        
        /**
         * loop main bam file 
         * Force add 1 at bit 3 of sam flags (only on target seq)
         */
        SAMRecordIterator itMain = samReaderInput.iterator();
  
        while(itMain.hasNext()){
            SAMRecord data = itMain.next();
           
            String seq = data.getReadString();
            MessageDigest md = MessageDigest.getInstance("MD5");
            md.update(seq.getBytes());
            seqMd5 = md.digest();
            String strSeqMd5 = DatatypeConverter.printHexBinary(seqMd5).toUpperCase();                

            if(targetSeqMd5Map.containsKey(strSeqMd5)){
                int numFlag = data.getFlags();
                if(numFlag == 3 || numFlag==2){
                    StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                    strFlag.append('1');
                    data.setFlags(Integer.parseInt(strFlag.toString(),2));
                }else if(numFlag==1 || numFlag==0){
                    StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                    strFlag.append("10");
                    data.setFlags(Integer.parseInt(strFlag.toString(),2));
                }else{
                    StringBuilder strFlag = new StringBuilder(Integer.toBinaryString(numFlag));
                    strFlag.setCharAt((strFlag.length()-1)-2,'1');
                    data.setFlags(Integer.parseInt(strFlag.toString(),2));
                }
                count++;
            }
            outputSam.addAlignment(data);
        }
        /*****************************************************/
        
        outputSam.close();
        samReaderInput.close();
        samReaderTarget.close();
        
        System.out.println("num read modifued = "+count);
    }
}
