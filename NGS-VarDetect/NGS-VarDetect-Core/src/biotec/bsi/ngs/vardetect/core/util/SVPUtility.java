/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.CombineReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.ThreadPool;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Writer;
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
    
    public static void fastIntersectCSVReport(String fileA, String fileB, String outputPath, String svType) throws FileNotFoundException, IOException{
        // This function will return 3 report file 1. unique SV report of fileA 2. unique SV report of fileB 3. intersersect report of file A and B
        // user must define type of sv ["DEL","TAN","CHI","INTER","INTRA"]
        
        File inputFileA = new File(fileA);
        File inputFileB = new File(fileB);
        String fileAName = inputFileA.getName().split("\\.csv")[0];
        String fileBName = inputFileB.getName().split("\\.csv")[0];
        String headerA = "";
        String headerB = "";
        
        boolean eof = false;
        
        Map<String,String> fileAMap = new LinkedHashMap();
        Map<String,String> fileBMap = new LinkedHashMap();
        Map<String,String> uniqueMapA = new LinkedHashMap();
        Map<String,String> uniqueMapB = new LinkedHashMap();
        Map<String,String> intersectMap = new LinkedHashMap();
        
        BufferedReader br = null;
        String line = "";
        String cvsSplitBy = ",";
        
        if(svType.equals("DEL")){
            
        }
        
        
        // Read fileA
        try {

            br = new BufferedReader(new FileReader(fileA));
            headerA = br.readLine();
            while ((line = br.readLine()) != null) {

                // use comma as separator
                String[] svInfo = line.split(cvsSplitBy);
                
                if(svType.equals("DEL")||svType.equals("TAN")||svType.equals("CHI")){
                    String key = svInfo[1]+","+svInfo[2]+","+svInfo[3]+","+svInfo[4]+","+svInfo[5]+","+svInfo[6];
                    fileAMap.put(key, line);
                }else if(svType.equals("INTER")||svType.equals("INTRA")){
                    String key = svInfo[1]+","+svInfo[2]+","+svInfo[3]+","+svInfo[4]+","+svInfo[5]+","+svInfo[6]+","+svInfo[13]+","+svInfo[14]+","+svInfo[15]+","+svInfo[16]+","+svInfo[17]+","+svInfo[18];
                    fileAMap.put(key, line);
                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        
        // Read fileB
        try {

            br = new BufferedReader(new FileReader(fileB));
            headerB = br.readLine();
            while ((line = br.readLine()) != null) {

                // use comma as separator
                String[] svInfo = line.split(cvsSplitBy);
                
                if(svType.equals("DEL")||svType.equals("TAN")||svType.equals("CHI")){
                    String key = svInfo[1]+","+svInfo[2]+","+svInfo[3]+","+svInfo[4]+","+svInfo[5]+","+svInfo[6];
                    fileBMap.put(key, line);
                }else if(svType.equals("INTER")||svType.equals("INTRA")){
                    String key = svInfo[1]+","+svInfo[2]+","+svInfo[3]+","+svInfo[4]+","+svInfo[5]+","+svInfo[6]+","+svInfo[13]+","+svInfo[14]+","+svInfo[15]+","+svInfo[16]+","+svInfo[17]+","+svInfo[18];
                    fileBMap.put(key, line);
                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        
        // extract unique Map and intersect map
        uniqueMapA = mapDifference(fileAMap,fileBMap);
        uniqueMapB = mapDifference(fileBMap,fileAMap);
        intersectMap = getMapIntersection(fileAMap,fileBMap);
        
        // WriteFile
        String uniqueFileA = outputPath + File.separator + fileAName + "_unique.csv";
        String uniqueFileB = outputPath + File.separator + fileBName + "_unique.csv";
        String intersectFile = outputPath + File.separator + fileAName + "_intersect_"+ fileBName + ".csv";
        
        Writer writer;
        // Write unique FileA
        File f = new File(uniqueFileA); //File object
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(uniqueFileA,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(uniqueFileA);
        }
        writer.write(headerA);
        writer.write("\n");        
        for(Map.Entry<String,String> entry : uniqueMapA.entrySet()){
            writer.write(entry.getValue());
            writer.write("\n");
        }
        writer.flush();
        writer.close();
        
        // Write unique FileB
        f = new File(uniqueFileB); //File object
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(uniqueFileB,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(uniqueFileB);
        }
        writer.write(headerB);
        writer.write("\n");        
        for(Map.Entry<String,String> entry : uniqueMapB.entrySet()){
            writer.write(entry.getValue());
            writer.write("\n");
        }
        writer.flush();
        writer.close();
        
        // Write intersect
        f = new File(intersectFile); //File object
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(intersectFile,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(intersectFile);
        }
        writer.write(headerA);
        writer.write("\n");        
        for(Map.Entry<String,String> entry : intersectMap.entrySet()){
            writer.write(entry.getValue());
            writer.write("\n");
        }
        writer.flush();
        writer.close();
    }
    
    public static <K, V> Map<K, V> mapDifference(Map<? extends K, ? extends V> left, Map<? extends K, ? extends V> right) {
        /**
         * Find different between 2 Map (left and right) (use key) 
         * function will return unique Map of left Map
         */
        Map<K, V> difference = new LinkedHashMap<>();
        difference.putAll(left);
        difference.putAll(right);
        difference.entrySet().removeAll(right.entrySet());
        return difference;
    }
    
    public static <K, V> Map<K, V> getMapIntersection(Map<? extends K, ? extends V> mapOne, Map<? extends K, ? extends V> mapTwo){
    Map intersection = new LinkedHashMap();
    for (Object key: mapOne.keySet()){
        if (mapTwo.containsKey(key))
           intersection.put(key, mapOne.get(key));
    }
    return intersection;
    }
}
