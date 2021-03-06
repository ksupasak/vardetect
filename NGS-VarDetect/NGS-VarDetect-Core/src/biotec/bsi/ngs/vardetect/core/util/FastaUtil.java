/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;



/**
 *
 * @author worawich
 */
public class FastaUtil {
    
    public static void createReferenceFromContig(String[] inputPath) throws IOException{
        
        /**
         * This function will combine set of reference contig file into one single reference file read contig from fasta file
         * Recreate it as a reference like chromosome file in human 
         * The amount of base limit for each chromosome is 28bit 
         * Also create index to trace back to each scaffold 
         */
        
        int mask28Bit = 268435455 ;
        int sliding = 60;
        int window = 60;
        
        File inputFile = new File(inputPath[0]);
        String path = inputFile.getParent();
        boolean smallContigStatus = false;
        
//        StringBuffer seq = new StringBuffer();
        Map<Long,String> index = new LinkedHashMap();
        ArrayList<String> listSeq = new ArrayList();
        ArrayList<Integer> chrList = new ArrayList();
        ArrayList<Integer> iniPosList = new ArrayList();
        ArrayList<Integer> lastPosList = new ArrayList();
        ArrayList<String> contigNameList = new ArrayList();
        ArrayList<String> refFileNameList = new ArrayList();
//        int iniPos = 0;
//        int lastPos = 0;
        int chrCount = 0;
//        boolean firstFlag = true;
//        String scaffoldName = null;
        
        File refFile = new File(path,"Ref_"+inputFile.getName()); //File object
        FileWriter writer;
        if(refFile.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(refFile,false);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(refFile);
        }


        for(int numF=0;numF<inputPath.length;numF++){
            chrCount = chrCount+1;
            inputFile = new File(inputPath[numF]);     
            String name = inputFile.getName().split("\\.")[0];

            StringBuffer seq = new StringBuffer();
            int iniPos = 0;
            int lastPos = 0;
            boolean firstFlag = true;
            String scaffoldName = null;
            /**
             * Read file
             */
            try(BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(inputFile)))){
                String sb = null;
                while ((sb = reader.readLine()) != null){

                    if(sb.charAt(0)=='>'){


                        if(firstFlag == false){
                            lastPos = seq.length()-1;           // index is start at 0 So, we have to minus 1 in order to get actual index

    //                        long indexCode = ((long)chrCount<<28|(long)iniPos)<<28|(long)lastPos;       // iniPos and lastPos is actual index. it can be used directly
    //                        
    //                        index.put(indexCode,scaffoldName);

                            chrList.add(chrCount);
                            iniPosList.add(iniPos);
                            lastPosList.add(lastPos);
                            contigNameList.add(scaffoldName);
                            refFileNameList.add(name);

                            for(int i=0;i<1000;i++){
                               /**
                                * append "N" 1000 character
                                */
                               seq.append("N");
                            }

                            iniPos = seq.length();              // index is start at 0, the next iniPos is actual index plus 1 but to get actual index we have to minus 1. So, both cancel out and we don't have to do any thing.

                        }

                        scaffoldName = sb.substring(1);
                        firstFlag = false;

                    }else{

                        if(sb.length()>18){

                            if(seq.length()+sb.trim().length() > mask28Bit){

                                /**
                                 * Reach maximum number of base in one seq
                                 * initiate new seq and also reset iniPos to zero and firstFlag to true
                                 * Reset or Update: 
                                 *  iniPos to zero 
                                 *  firstFlag to true
                                 *  chrCount plus one
                                 */

                                listSeq.add(seq.toString());
                                seq = new StringBuffer();
                                seq.append(sb.trim());
                                
                                
                                
                                chrCount++;
    //                            firstFlag = true;
                                iniPos = 0;                            
                            }else{
                                seq.append(sb.trim());
                                
                                
                                
                            }
                            smallContigStatus = false;
                        }else{
                            firstFlag = true;
                            smallContigStatus = true;
                        }
                    } 
                }

                /**
                 * for last contig
                 * But add only long contig not small contig (<18 bp)
                 * Can check from smallContigStatus
                 */
                
                if(smallContigStatus == false){
                    lastPos = seq.length()-1;           // index is start at 0 So, we have to minus 1 in order to get actual index

                    chrList.add(chrCount);
                    iniPosList.add(iniPos);
                    lastPosList.add(lastPos);
                    contigNameList.add(scaffoldName);
                    refFileNameList.add(name);
                    listSeq.add(seq.toString());
                }
                /*******************************/

            }catch(IOException e){

            }
        }
        
        
//        File refFile = new File(path,"ContigRef.fa"); //File object
//        FileWriter writer;
//        if(refFile.exists()){
////            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(refFile,false);
//        }else{
////            ps = new PrintStream(filename);
//            writer = new FileWriter(refFile);
//        }
        
        /**
         * Write reference file
         */
        
//        for(int i=0;i<listSeq.size();i++){
//            int numChr = i+1;
//            writer.write(">chr"+numChr+"\n");
//            String dummyString = listSeq.get(i);
//            writer.write(listSeq.get(i));           
//        }
//        writer.flush();
//        writer.close();
        for(int i=0;i<listSeq.size();i++){
            int numChr = i+1;
            writer.write(">chr"+numChr+"\n");
            String dummyString = listSeq.get(i);
            int n = (dummyString.length()-window)/sliding;       
            for(int numWin=0;numWin<n;numWin++){
                String s = dummyString.substring(numWin*sliding,numWin*sliding+window);
                writer.write(s+"\n");
            }
        }
        writer.flush();
        writer.close();
        
        
        
        
        
        
        File f = new File(path,"Ref_"+inputFile.getName().split("\\.")[0]+".index"); //File object
        FileWriter writerIdx;
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writerIdx = new FileWriter(f,false);
        }else{
//            ps = new PrintStream(filename);
            writerIdx = new FileWriter(f);
        }

        /**
         * write index file
         */

        for(int i=0;i<chrList.size();i++){
            writerIdx.write(String.format("%s,%s,chr%d,%d,%d", refFileNameList.get(i),contigNameList.get(i),chrList.get(i),iniPosList.get(i),lastPosList.get(i)));
            writerIdx.write("\n");
        }
        writerIdx.flush();
        writerIdx.close();

    }
    
    public static void filterSampleFile(String inputPath) throws IOException{
        int count = 0;
        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(inputPath);
        String name = null;
        boolean forceBreakFlag =  false;
    
        StringBuilder seq = new StringBuilder();
        
        String filename = path.getParent()+File.separator+path.getFileName().toString().split("\\.")[0]+"_filter."+path.getFileName().toString().split("\\.")[1];
        PrintStream ps;
        FileWriter writer;        
        /**
         * Check File existing
         */
        
        File f = new File(filename); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(filename,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(filename);
        }
        

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    name = null;
                }else{
                    if(line.charAt(0)=='>'){
                        count++;
                        if(seq.length()>0){
                            writer.write(">"+name+"\n");
                            writer.write(seq+"\n");
                            seq = new StringBuilder();
                        }
                        name = line.substring(1);

                    }else if(line.charAt(0)=='N'){
                        name = null;
                        seq = new StringBuilder();
                    }else{                                           
                        seq.append(line.toString()); 
                    }
                }    
            }                      
            /**
             * Add data for last seq of a file
             * in order to check it is last seq of file or last seq from force break
             * Can check from forceBreakFlag;
             */
            
            
        }
        
        // for last line
        if(seq.length()>0){
            writer.write(">"+name+"\n");
            writer.write(seq+"\n");   
        }

        writer.flush();
        writer.close();
        
    }
    
    public static void reIndexChrNameFastaFile(String inputPath) throws IOException{
        /**
         * This function will reIndex of chromosome name which has non numeric representation to numeric number
         * Save new file with name _reIndex
         * And also save the index file for mapping back to original chr representation
         */
        
        int count = 0;
        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(inputPath);
        String name = null;
        boolean forceBreakFlag =  false;
    
        StringBuilder seq = new StringBuilder();
        
        String filename = path.getParent()+File.separator+path.getFileName().toString().split("\\.")[0]+"_reIndex.fa";
        String filenameIndex = path.getParent()+File.separator+path.getFileName().toString().split("\\.")[0]+"_reIndex.index";
        PrintStream ps;
        FileWriter writer;  
        FileWriter writerIndex;
        /**
         * Check File existing
         */
        
        File f = new File(filename); //File object        
        if(f.exists()){
            // file is already exist no need to do it again
//            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(filename,true);
            return;
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(filename);
        }
        writerIndex = new FileWriter(filenameIndex);

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)=='>'){
                        count++;
                        writer.write(">chr"+count+"\n");
                        writerIndex.write(path.getFileName()+","+line.substring(1)+",chr"+count+"\n");
                    }else{                                           
                        writer.write(line+"\n"); 
                    }
                }    
            }
        }
        
        writer.flush();
        writer.close();
        writerIndex.flush();
        writerIndex.close();
        
    }
    
    public static void reIndexChrNameFastaFileV2(String inputPath) throws IOException{
        /**
         * This function will reIndex of chromosome name which has non numeric representation to numeric number
         * Save new file with name _reIndex
         * And also save the index file for mapping back to original chr representation
         * 
         * The version 2 function will consider the number (not finish)
         */
        
        int count = 0;
        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(inputPath);
        String name = null;
        boolean forceBreakFlag =  false;
    
        StringBuilder seq = new StringBuilder();
        
        String filename = path.getParent()+File.separator+path.getFileName().toString().split("\\.")[0]+"_reIndex.fa";
        String filenameIndex = path.getParent()+File.separator+path.getFileName().toString().split("\\.")[0]+"_reIndex.index";
        PrintStream ps;
        FileWriter writer;  
        FileWriter writerIndex;
        /**
         * Check File existing
         */
        
        File f = new File(filename); //File object        
        if(f.exists()){
            // file is already exist no need to do it again
//            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(filename,true);
            return;
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(filename);
        }
        writerIndex = new FileWriter(filenameIndex);

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)=='>'){
                        count++;
                        writer.write(">chr"+count+"\n");
                        writerIndex.write(path.getFileName()+","+line.substring(1)+",chr"+count+"\n");
                    }else{                                           
                        writer.write(line+"\n"); 
                    }
                }    
            }
        }
        
        writer.flush();
        writer.close();
        writerIndex.flush();
        writerIndex.close();
        
    }
    
    public static void createSampleFromAlignResult(String alignResultFile, String sampleFile, char option) throws IOException{
        /**
         * This function will create Sample from alignment Result.
         * Option:
         *      'c' = It will cut out the front part of original Sample at specific index.Each sample have different cutting index. The cutting index is the first index that match on reference
         *            the cutting can get from alignment Result. You must use alignResult_Sorted.txt as a input.
         *      'r' = got complete sequence (not cut)
         * 
         */
        String cutSampleFilename = null;
        Map<String,String> cutSampleList = new LinkedHashMap();
        int count = 0;
        Charset charset = Charset.forName("US-ASCII");
        Path pathResult = Paths.get(alignResultFile);
        Path pathSample = Paths.get(sampleFile);
        String name = null;
        boolean forceBreakFlag =  false;
        
        StringBuilder seq = new StringBuilder();
        
        if(option == 'c'){
            cutSampleFilename = pathSample.getParent()+File.separator+pathSample.getFileName().toString().split("\\.")[0]+"_frontCutSample."+pathSample.getFileName().toString().split("\\.")[1];
        }else if(option == 'r'){
            cutSampleFilename = pathSample.getParent()+File.separator+pathSample.getFileName().toString().split("\\.")[0]+"_Sample."+pathSample.getFileName().toString().split("\\.")[1];
        }
        
        
        PrintStream ps;
        FileWriter writer;  
        
        /**
         * Check File existing
         */
        
        File f = new File(cutSampleFilename); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(cutSampleFilename,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(cutSampleFilename);
        }
        
        /**
         * Read alignment result
         * set to Arraylist
         */
        ArrayList<String> resultList = new ArrayList();
        
        try (BufferedReader reader = Files.newBufferedReader(pathResult, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {                
                resultList.add(line);                
            }
        }
        
        /**
         * Read sample file
         * set to Map
         */
        Map<String,String> sampleList = new LinkedHashMap();
        try (BufferedReader reader = Files.newBufferedReader(pathSample, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {                
                if(line.charAt(0)=='>'){
                    if(seq.length()>0){
                        sampleList.put(name, seq.toString());
                        seq = new StringBuilder();
                    }
                    name = line.substring(1);
                }else{                                           
                    seq.append(line.toString()); 
                }                
            }
        }
        
        if(option == 'c'){
            /**
            * Cut sample
            * write to file
            */
            
            Map<String,Integer> sampleCutPoint = findSampleCutPoint(resultList);
        
        
            for (Map.Entry<String, Integer> entry : sampleCutPoint.entrySet()){
                String readName = entry.getKey();
                int cutPoint = entry.getValue();

                String sequence = sampleList.get(readName);
                String cutSequence = sequence.substring(0, cutPoint);
                if(cutPoint != 0){
                    writer.write(">"+readName+"\n");
                    writer.write(cutSequence+"\n");
                }               
    //            cutSampleList.put(name, cutSequence);
            }
            
        }else if(option == 'r'){
            /**
             * Not cut sample
             * write to file
             */
            
            Map<String,Integer> sampleCutPoint = findSampleCutPoint(resultList);
        
        
            for (Map.Entry<String, Integer> entry : sampleCutPoint.entrySet()){
                String readName = entry.getKey();

                String sequence = sampleList.get(readName);

                writer.write(">"+readName+"\n");
                writer.write(sequence+"\n");
    //            cutSampleList.put(name, cutSequence);
            }
        }
        
        writer.flush();
        writer.close();    
    }
    
    public static void createSampleFromNonVariantCoverageReport(String alignResultFile, String sampleFile) throws IOException{
        /**
         * This function will create Sample from Coverage report.
         * Separately save to file per group
         * Save to fasta format (.fa)
         */
        String cutSampleFilename = null;
       
        Charset charset = Charset.forName("US-ASCII");
        Path pathResult = Paths.get(alignResultFile);
        Path pathSample = Paths.get(sampleFile);
        String name = null;
        StringBuilder seq = new StringBuilder();
        
        
        
        
        
        /**
         * Read coverage report
         * set to groupMap
         */
        int groupNumber = 0;
        Map<Integer,ArrayList<String>> groupMap = new LinkedHashMap(); // Map contain key  is group number and value is list of read name in that group.
        
        try (BufferedReader reader = Files.newBufferedReader(pathResult, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {
                String[] splitData = line.split("\\s+");
                if(splitData[0].equals("Group")){
                    groupNumber = Integer.parseInt(splitData[2]);
                }else{
                    splitData = line.split(",");
                    String readName = splitData[9];
                    
                    if(groupMap.containsKey(groupNumber)){
                        ArrayList<String> readNameList = groupMap.get(groupNumber);
                        readNameList.add(readName);
                        groupMap.put(groupNumber, readNameList);
                    }else{
                        ArrayList<String> readNameList = new ArrayList();
                        readNameList.add(readName);
                        groupMap.put(groupNumber, readNameList);
                    }
                }                 
            }
        }
        /***************************************************/
        
        /**
         * Read sample file
         * set to Map
         */
        Map<String,String> sampleList = new LinkedHashMap();  // key is name and value is sequence
        try (BufferedReader reader = Files.newBufferedReader(pathSample, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {                
                if(line.charAt(0)=='>'){
                    if(seq.length()>0){
                        sampleList.put(name, seq.toString());
                        seq = new StringBuilder();
                    }
                    name = line.substring(1);
                }else{                                           
                    seq.append(line.toString()); 
                }                
            }
        }
        /***************************************************/
        
        /**
         * extract and write to file
         */
        
        for (Map.Entry<Integer, ArrayList<String>> entry : groupMap.entrySet()){
            groupNumber = entry.getKey();
            ArrayList<String> readNameList = entry.getValue();
            
            cutSampleFilename = pathSample.getParent()+File.separator+pathSample.getFileName().toString().split("\\.")[0]+"Group_"+groupNumber+"_Sample.fa";
            PrintStream ps;
            FileWriter writer;  

            /**
             * Check File existing
             */

            File f = new File(cutSampleFilename); //File object        
            if(f.exists()){
    //            ps = new PrintStream(new FileOutputStream(filename,true));
                writer = new FileWriter(cutSampleFilename,true);
            }else{
    //            ps = new PrintStream(filename);
                writer = new FileWriter(cutSampleFilename);
            }
            
            
            for(int i=0;i<readNameList.size();i++){
                String readName = readNameList.get(i);
                String sequence = sampleList.get(readName);
                
                writer.write(">"+readName+"\n");
                writer.write(sequence+"\n");
            }
            
            writer.flush();
            writer.close();
        }
        /*******************************/
   
    }
    
    public static void createSampleFromVariantCoverageReport(String coverageResultFile, String sampleFile) throws IOException{
        /**
         * This function will create Sample from Coverage report.
         * Separately save to file per group
         * Support fasta and fastq sample file
         * Save to fasta format (.fa)
         */
        String cutSampleFilename = null;
       
        Charset charset = Charset.forName("US-ASCII");
        Path pathResult = Paths.get(coverageResultFile);
        Path pathSample = Paths.get(sampleFile);
        String name = null;
        StringBuilder seq = new StringBuilder();
        
        String[] strdummy = sampleFile.split("\\.");
        String sampleFileType = strdummy[strdummy.length-1];
        
        /**
         * Read coverage report
         * set to groupMap
         */
        int groupNumber = 0;
        Map<Integer,ArrayList<String>> groupMap = new LinkedHashMap(); // Map contain key  is group number and value is list of read name in that group.
        
        try (BufferedReader reader = Files.newBufferedReader(pathResult, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {
                String[] splitData = line.split("\\s+");
                if(splitData[0].equals("Group")){
                    groupNumber = Integer.parseInt(splitData[1]);
                }else if(!splitData[0].equals("Variation")){
                    splitData = line.split(",");
                    String readName = splitData[9];
                    if(readName.equals("null")){
                        // in one tail at back part. cause front part has null read name
                        readName = splitData[21];   // index 21 is index of backpart read name
                    }
                    
                    if(groupMap.containsKey(groupNumber)){
                        ArrayList<String> readNameList = groupMap.get(groupNumber);
                        readNameList.add(readName);
                        groupMap.put(groupNumber, readNameList);
                    }else{
                        ArrayList<String> readNameList = new ArrayList();
                        readNameList.add(readName);
                        groupMap.put(groupNumber, readNameList);
                    }
                }                 
            }
        }
        /***************************************************/
        
        /**
         * Read sample file
         * set to Map
         */
        Map<String,String> sampleList = new LinkedHashMap();  // key is name and value is sequence
        if(sampleFileType.endsWith("fa")||sampleFileType.endsWith("fasta")){            
            try (BufferedReader reader = Files.newBufferedReader(pathSample, charset)) {
                String line = null;                   
                while ((line = reader.readLine()) != null) {                
                    if(line.charAt(0)=='>'){
                        if(seq.length()>0){
                            sampleList.put(name, seq.toString());
                            seq = new StringBuilder();
                        }
                        name = line.substring(1);
                    }else{                                           
                        seq.append(line.toString()); 
                    }                
                }
            }
        }else if(sampleFileType.endsWith("fastq")){            
            try (BufferedReader reader = Files.newBufferedReader(pathSample, charset)) {
                String line = null;
                byte counter = 0;
                while ((line = reader.readLine()) != null) {
                    ++counter;
                    if(line.isEmpty()){
                        name = null;
                    }else{
                        if(counter==1){
                            if(seq.length()>0){                                
                                sampleList.put(name, seq.toString());
                                seq = new StringBuilder();
                            }
                            name = line.substring(1);
                            if(name.equals("ERR846998.3062537 HS25_15277:3:1301:12738:55130#59/2")){
                                System.out.println();
                            }
                        }else if(counter==2){                                               
                            seq.append(line.toString()); 
                        }
                    }
                    if(counter == 4){
                        counter=0;
                    }
                }
            }
        }
        
        /***************************************************/
        
        /**
         * extract and write to file
         */
        ArrayList<String> usedReadName = new ArrayList();
        cutSampleFilename = pathSample.getParent()+File.separator+pathResult.getFileName().toString().split("\\.")[0]+"_Sample.fa";
        PrintStream ps;
        FileWriter writer;  

        /**
         * Check File existing
         */

        File f = new File(cutSampleFilename); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(cutSampleFilename,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(cutSampleFilename);
        }
        
        for (Map.Entry<Integer, ArrayList<String>> entry : groupMap.entrySet()){
            groupNumber = entry.getKey();
            ArrayList<String> readNameList = entry.getValue();
            
            
            
            
            for(int i=0;i<readNameList.size();i++){
                String readName = readNameList.get(i);
                
                if(!usedReadName.contains(readName)){
                    String sequence = sampleList.get(readName);
                
                    writer.write(">"+readName+"\n");
                    writer.write(sequence+"\n");
                } 
            } 
        }
        writer.flush();
        writer.close();
        /*******************************/
   
    }
    
    public static Map<String,Integer> findSampleCutPoint(ArrayList<String> inResultList){
        String oldReadName = null;
        Map<String,Integer> sampleCutPoint = new LinkedHashMap();
        ArrayList<Integer> iniIndexList = new ArrayList();
        boolean firstTimeFlag = true;
        
        for(int i=0;i<inResultList.size();i++){
            String dataGet = inResultList.get(i);
            // continue this point find iniindex to use ass cut point
            /***    Extract data    ****/
            String[] data = dataGet.split(",");           
            int iniIdx = Integer.parseInt(data[8]);
            String readName = data[9];
            /******************************/

            if(readName.equals(oldReadName)){
                /* Same set of read (continue add data) */
                
                iniIndexList.add(iniIdx);
                
            }else{
                /**
                 * Found new set of Read 
                 * 1. find min iniIndex
                 * 2. put to map
                 * 3. create new iniIndexList       
                 */
                if(firstTimeFlag == false){
                    int min = iniIndexList.get(0);
                    for(Integer num: iniIndexList){
                        if(num < min) min = num;
                    }
                    
                    sampleCutPoint.put(oldReadName, min);
                    
                    iniIndexList = new ArrayList();
                }  
            }
            
            iniIndexList.add(iniIdx);
            oldReadName = readName;
            firstTimeFlag = false; 
        }
        
        return sampleCutPoint;
    }
    
    public static void separateContigToFile(String filename) throws IOException{
        Path pathFile = Paths.get(filename);
        
        String savePath = pathFile.getParent()+File.separator;
        Charset charset = Charset.forName("US-ASCII");
        FileWriter writer = null;
        boolean firstTimeFlag = true;
               
        try (BufferedReader reader = Files.newBufferedReader(pathFile, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {                
                if(line.charAt(0)== '>'){
                    if(firstTimeFlag == false){
                        writer.flush();
                        writer.close();
                    }
                    firstTimeFlag = false;
                    String saveFilename = line.substring(1);
                    writer = new FileWriter(savePath+"contig_"+saveFilename+".fa");
                    writer.write(line+"\n");
                }else{
                    writer.write(line+"\n");
                }           
            }
        }
        writer.flush();
        writer.close();
    }
    
    public static void createTargetEventFromRead(String inputPath) throws IOException{
        
        ///**** Not finish (ยังไม่ได้ทำอะไรแค่ก๊อปโครงมาจะเอาส่วนอ่านไฟล์มาใช้เฉยๆ)
        
        /**
         * This function will reIndex of chromosome name which has non numeric representation to numeric number
         * Save new file with name _reIndex
         * And also save the index file for mapping back to original chr representation
         */
        
        int count = 0;
        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(inputPath);
        String name = null;
        boolean forceBreakFlag =  false;
    
        StringBuilder seq = new StringBuilder();
        
        String filename = path.getParent()+File.separator+path.getFileName().toString().split("\\.")[0]+"_reIndex.fa";
        String filenameIndex = path.getParent()+File.separator+path.getFileName().toString().split("\\.")[0]+"_reIndex.index";
        PrintStream ps;
        FileWriter writer;  
        FileWriter writerIndex;
        /**
         * Check File existing
         */
        
        File f = new File(filename); //File object        
        if(f.exists()){
            // file is already exist no need to do it again
//            ps = new PrintStream(new FileOutputStream(filename,true));
//            writer = new FileWriter(filename,true);
            return;
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(filename);
        }
        writerIndex = new FileWriter(filenameIndex);

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;                   
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)=='>'){
                        count++;
                        writer.write(">chr"+count+"\n");
                        writerIndex.write(path.getFileName()+","+line.substring(1)+",chr"+count+"\n");
                    }else{                                           
                        writer.write(line+"\n"); 
                    }
                }    
            }
        }
        
        writer.flush();
        writer.close();
        writerIndex.flush();
        writerIndex.close();
        
    }
    
    public static void createIntesectFastaFile(String fa1, String fa2) throws IOException{
        /**
         * This function will find the intersect read sequence between two fasta file
         * The smallest file will be pick as main file store sequence in Map
         * walk through another file and chack intersect
         * export intersect read sequence in fasta
         */
        
        File fasta1 = new File(fa1);
        File fasta2 = new File(fa2);
        String exportFile = "";
        String mainFile = "";
        String minorFile = "";
        
        Map<String,String> mainSeqMap = new LinkedHashMap();
        Map<String,String> intersectSeqMap = new LinkedHashMap();
        
        if(fasta1.length() > fasta2.length()){
            mainFile = fa2;
            minorFile = fa1;
            
            exportFile = fasta2.getParent()+ File.separator + fasta2.getName().split("\\.")[0] + "_Intersect_" + fasta1.getName();
        }else{
            mainFile = fa1;
            minorFile = fa2;
            
            exportFile = fasta1.getParent() + File.separator + fasta1.getName().split("\\.")[0] + "_Intersect_" + fasta2.getName();
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
                    if(line.charAt(0)=='>'){                        
                        name = line.substring(1);
                    }else{                                           
                        mainSeqMap.put(name,line);
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
            boolean intersectFlag = false;
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)=='>'){                        
                        name = line.substring(1);
                        if(mainSeqMap.containsKey(name)){
                            intersectFlag = true;
                        }else{
                            intersectFlag = false;
                        }
                    }else{
                        if(intersectFlag == true){
                            intersectSeqMap.put(name, line);
                        }
                    }
                }    
            }
        }        
        /************************************************************/
        
        /**
         * Export intersect read in fasta file format
         */
        
        FileWriter writer;
        writer = new FileWriter(exportFile);

        for(Map.Entry<String,String> entry: intersectSeqMap.entrySet()){ 
            writer.write(">"+entry.getKey());
            writer.write("\n");
            writer.write(entry.getValue());
            writer.write("\n");
        }
        /**************************************************************/
    }
    
    public static void createUnIntesectFastaFile(String fa1, String fa2) throws IOException{
        /**
         * This function will find the intersect read sequence between two fasta file
         * The smallest file will be pick as main file store sequence in Map
         * walk through another file and chack intersect
         * export intersect read sequence in fasta
         */
        
        File fasta1 = new File(fa1);
        File fasta2 = new File(fa2);
        String exportFile = "";
        String mainFile = "";
        String minorFile = "";
        
        Map<String,String> mainSeqMap = new LinkedHashMap();
        Map<String,String> unIntersectSeqMap = new LinkedHashMap();
        
        if(fasta1.length() > fasta2.length()){
            mainFile = fa2;
            minorFile = fa1;
            
            exportFile = fasta2.getParent()+ File.separator + fasta2.getName().split("\\.")[0] + "_unIntersect_" + fasta1.getName();
        }else{
            mainFile = fa1;
            minorFile = fa2;
            
            exportFile = fasta1.getParent() + File.separator + fasta1.getName().split("\\.")[0] + "_unIntersect_" + fasta2.getName();
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
                    if(line.charAt(0)=='>'){                        
                        name = line.substring(1);
                    }else{                                           
                        mainSeqMap.put(name,line);
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
            boolean unIntersectFlag = false;
            while ((line = reader.readLine()) != null) {
                if(line.isEmpty()){
                    
                }else{
                    if(line.charAt(0)=='>'){                        
                        name = line.substring(1);
                        if(mainSeqMap.containsKey(name)){
                            unIntersectFlag = false;
                        }else{
                            unIntersectFlag = true;
                        }
                    }else{
                        if(unIntersectFlag == true){
                            unIntersectSeqMap.put(name, line);
                        }
                    }
                }    
            }
        }        
        /************************************************************/
        
        /**
         * Export intersect read in fasta file format
         */
        
        FileWriter writer;
        writer = new FileWriter(exportFile);

        for(Map.Entry<String,String> entry: unIntersectSeqMap.entrySet()){ 
            writer.write(">"+entry.getKey());
            writer.write("\n");
            writer.write(entry.getValue());
            writer.write("\n");
        }
        /**************************************************************/
    }

}
