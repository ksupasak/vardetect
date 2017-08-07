/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.LASTFusionPattern;
import biotec.bsi.ngs.vardetect.core.LASTResult;
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
public class LASTUtils {
    
    public static void countSampleFromLASTResult(String inLASTtResultFile) throws IOException{
        /**
         * This function will read the result file from LAST
         * Then count the number of unique read and count the number of location match of each read 
         * Export to file _LASTCountSample.txt
         */
        
        Charset charset = Charset.forName("US-ASCII");
        //String[] ddSS = filename.split(".");
        String saveFileName = inLASTtResultFile.split("\\.")[0] + "_LASTCountSample.txt";
        Path path = Paths.get(inLASTtResultFile);
        Map<String,Integer> mapCountSample = new LinkedHashMap();
        int dummy = 1;
        String sampleName = "";
        String data = "";
        ArrayList<String> inData = new ArrayList();    
        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            
            while ((data = reader.readLine()) != null) {
                if(data.isEmpty() == false){ 
                    if(data.charAt(0) == 's' && dummy == 1){
                        dummy++;
                    }else if(data.charAt(0) == 's' && dummy == 2) {
                        dummy = 1;
                        String[] splitData = data.split("\\s+");
                        sampleName = splitData[1];

                        if(mapCountSample.containsKey(sampleName)){
                            int amount = mapCountSample.get(sampleName);
                            mapCountSample.put(sampleName, ++amount);
                        }else{
                            mapCountSample.put(sampleName, 1);
                        }
                    }                
                }
            }                
        }
        
        FileWriter writer;        
        /**
         * Check File existing
         */
        ArrayList<Integer> countRead = new ArrayList();
        File f = new File(saveFileName); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(saveFileName,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(saveFileName);
        }
        
        writer.write("Amount of align read : "+mapCountSample.size()+"\n");
        writer.write("\nList of sample and number of location match\n");
        for(Map.Entry<String,Integer> entry : mapCountSample.entrySet()){
            String name = entry.getKey();
            int number = entry.getValue();
            
            writer.write("Read Name : "+name+"\tNumber of location match : "+number+"\n");         
        }
        writer.flush();
        writer.close();
    }
    
    public static void countIntersectFromLASTResult(String lastResult1, String lastResult2) throws IOException{
        /**
         * This function will loop over two file lastResult1 and lastResult2
         * this function will use lastResult1 as a main file and let lastResult2 to compare with it
         * It mean it will read and store data of lastResult1 in memory but not for lastResult2
         * the program will read through lastResult2 and compare with the data in memory simultaneously
         * 
         * Before start this function will check the sample filename of two input. It will start find intersect when there is the same sample filename
         * The file name is any character that came before "." but not include path (if it has several "." it will pick the first part before first ".")
         */
        String[] splitData1 = lastResult1.split("/");
        String[] splitData2 = lastResult2.split("/");
        String sampleFileName1 = splitData1[splitData1.length-1].split("\\.")[0];
        String sampleFileName2 = splitData2[splitData2.length-1].split("\\.")[0];
        
        if(sampleFileName1.equals(sampleFileName2)){            
        
            ArrayList<LASTResult> listLastResult1 = readLASTResult(lastResult1);
            Map<String,ArrayList<LASTResult>> readNameCheckMap = new LinkedHashMap();
            Map<String,ArrayList<LASTResult>> readNameIntersectMap = new LinkedHashMap();
            /**
             * loop to create map for checking contain readName in lastResult1
             */
            for(LASTResult object : listLastResult1){
                String sName = object.getSampleName();

                if(readNameCheckMap.containsKey(sName)){
                    ArrayList<LASTResult> dummyList = readNameCheckMap.get(sName);
                    dummyList.add(object);
                    readNameCheckMap.put(sName, dummyList);      
                }else{
                    ArrayList<LASTResult> dummyList = new ArrayList();
                    dummyList.add(object);
                    readNameCheckMap.put(sName, dummyList);
                }           
            }

            /**
             * read file lastResult2 and check contain simultaneously
             */
            Charset charset = Charset.forName("US-ASCII");

            int score=0;
            String EG2="";
            String E="";

            String refChr="";
            long refPos=0;
            int refNumBase=0;
            String refStrand="";
            int refSize=0;
            String refSequenceMatch="";
            long refSequenceCode=0;

            String sampleName="";
            long samplePos=0;
            int sampleNumBase=0;
            String sampleStrand="";
            int sampleSize=0;
            String sampleSequenceMatch="";
            long sampleSequenceCode=0;

            int dummy = 1;
            String data = "";

            Path path = Paths.get(lastResult2);
            try (BufferedReader reader = Files.newBufferedReader(path, charset)) {

                while ((data = reader.readLine()) != null) {
                    if(data.isEmpty() == false){ 
                        if(data.charAt(0) == 's' && dummy == 1){
                            dummy++;
                            String[] splitData = data.split("\\s+");
                            refChr = splitData[1];
                            refPos = Long.parseLong(splitData[2]);
                            refNumBase = Integer.parseInt(splitData[3]);
                            refStrand = splitData[4];
                            refSize = Integer.parseInt(splitData[5]);
                            refSequenceMatch = splitData[6];
                            refSequenceCode = SequenceUtil.encodeMer(refSequenceMatch, 18);

                        }else if(data.charAt(0) == 's' && dummy == 2) {
                            dummy = 1;

                            String[] splitData = data.split("\\s+");
                            sampleName = splitData[1];                        
                            samplePos = Long.parseLong(splitData[2]);
                            sampleNumBase = Integer.parseInt(splitData[3]);
                            sampleStrand = splitData[4];
                            sampleSize = Integer.parseInt(splitData[5]);
                            sampleSequenceMatch = splitData[6];
                            sampleSequenceCode = SequenceUtil.encodeMer(sampleSequenceMatch, 18); 

                            if(readNameCheckMap.containsKey(sampleName)){
                                LASTResult lastRes = new LASTResult();
                                lastRes.addRefData(refChr,refPos,refNumBase,refStrand,refSize,refSequenceMatch,refSequenceCode);
                                lastRes.addSampleData(sampleName,samplePos,sampleNumBase,sampleStrand,sampleSize,sampleSequenceMatch,sampleSequenceCode);
                                lastRes.addScoreData(score,EG2,E);

                                if(readNameIntersectMap.containsKey(sampleName)){
                                    ArrayList<LASTResult> dummy2 = readNameIntersectMap.get(sampleName);
                                    dummy2.add(lastRes);
                                    readNameIntersectMap.put(sampleName, dummy2);
                                }else{
                                    ArrayList<LASTResult> dummy2 = readNameCheckMap.get(sampleName);
                                    dummy2.add(lastRes);
                                    readNameIntersectMap.put(sampleName, dummy2);
                                }   
                            }
                        }else if(data.charAt(0)=='a'){
                            String[] splitData = data.split("\\s+");
                            score = Integer.parseInt(splitData[1].split("=")[1]);
                            EG2 = splitData[2].split("=")[1];
                            E = splitData[3].split("=")[1];
                        }                
                    }
                }                
            }
 
        
        
            /**
             * write result from readNameCheckMap
             */
            String[] dummyfilename2 = lastResult2.split("/");
            String sfilename2 = dummyfilename2[dummyfilename2.length-1].split("_")[1];
            String saveFileName = lastResult1.split("\\.")[0]+"_"+lastResult1.split("\\.")[1]+"_"+sfilename2+"_LASTIntesectResult.txt";
            FileWriter writer;        
            /**
             * Check File existing
             */
            ArrayList<Integer> countRead = new ArrayList();
            File f = new File(saveFileName); //File object        
            if(f.exists()){
    //            ps = new PrintStream(new FileOutputStream(filename,true));
                writer = new FileWriter(saveFileName,true);
            }else{
    //            ps = new PrintStream(filename);
                writer = new FileWriter(saveFileName);
            }

            writer.write("Amount of align read : "+readNameIntersectMap.size()+"\n");
            writer.write("\nList of LAST result group\n");
            int count = 0;
            for(Map.Entry<String,ArrayList<LASTResult>> entry : readNameIntersectMap.entrySet()){

                String samName = entry.getKey();
                ArrayList<LASTResult> listLASTResult = entry.getValue();
                writer.write("\nGroup:"+(++count)+"\tSequenceName:"+samName+"\tPattern:"+listLASTResult.size()+"\n");

                for(LASTResult LASTResultObject : listLASTResult){
                    writer.write(LASTResultObject.toString());
                }
            }
            writer.flush();
            writer.close();
        }
    }
    
    public static void analyzeFusionFromIntersecResult(String intersectFile, int threshold) throws IOException{
        /**
         * Read file and separate two database result and find every possible fusion between two database.
         * Group the same fusion pattern together then Filter low coverage (lower than threshold) group out and write to file (.txt)
         */
        
        Map<Integer,ArrayList<Integer>> lastResGroupMap = new LinkedHashMap();  // use to store lastresult object sperately by group. Key is group || Value is list of LAST result object ID 
        Map<String,Integer> chrIDMap = new LinkedHashMap();                     // use to store chr ID  Key is string of chr name || Value is chr ID
//        Map<Long,ArrayList<Integer>> fusionGroupMap = new LinkedHashMap();       // use to store group of fusion pattern  Key is fusionGroupID || Value is list of fusion pattern object ID
        Map<Long,Map<Long,ArrayList<LASTFusionPattern>>> fusionGroupMap = new LinkedHashMap();  // use to store group of fusion pattern. you have to pass to key to get to the array of fusion pattern first key is frontcode [chr|position] second is backCode [chr|position] last value is arrayList of object LASTFusionPattern. 
        ArrayList<LASTResult> listLASTResult = new ArrayList();                 // contaion LASTResult object. can use object ID to access data
        
        // continue feed data into map separate by group
        
        Charset charset = Charset.forName("US-ASCII");
        
       
        
        int score=0;
        String EG2="";
        String E="";
        
        String refChr="";
        long refPos=0;
        int refNumBase=0;
        String refStrand="";
        int refSize=0;
        String refSequenceMatch="";
        long refSequenceCode=0;
        
        String sampleName="";
        long samplePos=0;
        int sampleNumBase=0;
        String sampleStrand="";
        int sampleSize=0;
        String sampleSequenceMatch="";
        long sampleSequenceCode=0;
        
        int groupNumber = 0;
        int numPattern = 0;
        //String[] ddSS = filename.split(".");
//        String saveFileName = inLASTtResultFile.split("\\.")[0] + "_LASTCountSample.txt";
        Path path = Paths.get(intersectFile);
        Map<String,Integer> mapCountSample = new LinkedHashMap();
        int dummy = 1;
        int count = 0;
        int idCount = 0;
        String data = "";
        ArrayList<String> inData = new ArrayList();    
        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
           
            while ((data = reader.readLine()) != null) {
                if(data.isEmpty()==false){
                    if(data.charAt(0) == 'G'){
                        String[] dummyData = data.split("\\s+");
                        groupNumber = Integer.parseInt(dummyData[0].split(":")[1]);
                        numPattern = Integer.parseInt(dummyData[2].split(":")[1]);
//                        count = 0;
                    }else if(data.charAt(0) == 's' && dummy == 1){
                        dummy++;
                        String[] splitData = data.split("\\s+");
                        refChr = splitData[1];
                        refPos = Long.parseLong(splitData[2]);
                        refNumBase = Integer.parseInt(splitData[3]);
                        refStrand = splitData[4];
                        refSize = Integer.parseInt(splitData[5]);
                        refSequenceMatch = splitData[6];
                        refSequenceCode = SequenceUtil.encodeMer(refSequenceMatch, 18);
                    }else if(data.charAt(0) == 's' && dummy == 2) {
                        dummy = 1;

                        String[] splitData = data.split("\\s+");
                        sampleName = splitData[1];                        
                        samplePos = Long.parseLong(splitData[2]);
                        sampleNumBase = Integer.parseInt(splitData[3]);
                        sampleStrand = splitData[4];
                        sampleSize = Integer.parseInt(splitData[5]);
                        sampleSequenceMatch = splitData[6];
                        sampleSequenceCode = SequenceUtil.encodeMer(sampleSequenceMatch, 18); 

                        if(mapCountSample.containsKey(sampleName)){
                            int amount = mapCountSample.get(sampleName);
                            mapCountSample.put(sampleName, ++amount);
                        }else{
                            mapCountSample.put(sampleName, 1);
                        }

                        LASTResult lastRes = new LASTResult();
                        lastRes.addRefData(refChr,refPos,refNumBase,refStrand,refSize,refSequenceMatch,refSequenceCode);
                        lastRes.addSampleData(sampleName,samplePos,sampleNumBase,sampleStrand,sampleSize,sampleSequenceMatch,sampleSequenceCode);
                        lastRes.addScoreData(score,EG2,E);
                        listLASTResult.add(lastRes);

                        // Create chrID Map
                        if(!chrIDMap.containsKey(refChr)){
                            chrIDMap.put(refChr, ++idCount);
                        }
                        // add result data to lastResGroupMap (store data ID in map for later use)
                        if(lastResGroupMap.containsKey(groupNumber)){
                            ArrayList<Integer> listLASTResultID = lastResGroupMap.get(groupNumber);
                            listLASTResultID.add(count);
                            lastResGroupMap.put(groupNumber, listLASTResultID);  
                        }else{
                            ArrayList<Integer> listLASTResultID = new ArrayList();
                            listLASTResultID.add(count);
                            lastResGroupMap.put(groupNumber, listLASTResultID);
                        }
                        count++;
                    }else if(data.charAt(0)=='a'){
                        String[] splitData = data.split("\\s+");
                        score = Integer.parseInt(splitData[1].split("=")[1]);
                        EG2 = splitData[2].split("=")[1];
                        E = splitData[3].split("=")[1];
                    }
                }
            }
        }
        
        
        /**
         * Find possible fusion pattern
         */
        long frontCode = 0;
        long backCode = 0;
        for(Map.Entry<Integer,ArrayList<Integer>> entry : lastResGroupMap.entrySet()){
            String chrName = "";
            int numGroup = entry.getKey();
            ArrayList<Integer> listID = entry.getValue();
            for(int i=0;i < listID.size();i++){
                int mainID = listID.get(i);
                LASTResult mainLASTResult = listLASTResult.get(mainID);
                
                long mainIniIndex = mainLASTResult.getSamplePos();
                int mainNumBase = mainLASTResult.getSampleNumBase();
                long mainEndIndex = (mainIniIndex + mainNumBase)-1;
                long mainIniMatchPos = mainLASTResult.getRefPos();
                long mainEndMatchPos = (mainIniMatchPos + mainNumBase)-1;
                String mainChr = mainLASTResult.getRefChr();
                long mainChrID = chrIDMap.get(mainChr);
                
                for(int j=i+1;j<listID.size();j++){
//                    if(i==2){
//                        System.out.println();
//                    }
                    
                    int subID = listID.get(j);
                    LASTResult subLASTResult = listLASTResult.get(subID);
                    String subChr = subLASTResult.getRefChr();
                    if(subChr.equals(mainChr)){
                        continue;
                    }
                    long subIniIndex = subLASTResult.getSamplePos();
                    int subNumBase = subLASTResult.getSampleNumBase();
                    long subEndIndex = (subIniIndex + subNumBase)-1;
                    long subIniMatchPos = subLASTResult.getRefPos();
                    long subEndMatchPos = (subIniMatchPos + subNumBase)-1;
                    long subChrID = chrIDMap.get(subChr);
                    
                    LASTFusionPattern fusionPattern = new LASTFusionPattern();
                    
                    if(mainIniIndex < subIniIndex){                        
                        long overlapBase = (mainEndIndex-subIniIndex)+1;
                        frontCode = (mainChrID<<28)+mainEndMatchPos;
                        backCode = (subChrID<<28)+subIniMatchPos;                        
                        fusionPattern.setFrontID(mainID);
                        fusionPattern.setBackID(subID);
                        fusionPattern.setOverlapBase(overlapBase);
//                        if(mainID==subID){
//                            System.out.println();
//                        }
                    }else if(subIniIndex < mainIniIndex){
                        long overlapBase = (subEndIndex-mainIniIndex)+1;
                        frontCode = (subChrID<<28)+subEndMatchPos;
                        backCode = (mainChrID<<28)+mainIniMatchPos;
//                        if(mainID==subID){
//                            System.out.println();
//                        }
                        fusionPattern.setFrontID(subID);
                        fusionPattern.setBackID(mainID);
                        fusionPattern.setOverlapBase(overlapBase);
                    }else if(mainIniIndex == subIniIndex){
//                        // has same ini point
//                        if(subEndIndex>mainEndIndex){
//                            // sub is back part
//                            long overlapBase = (mainIniIndex - mainEndIndex);
//                            frontCode = (mainChrID<<28)+mainEndMatchPos;
//                            backCode = (subChrID<<28)+subIniMatchPos;
//                            fusionPattern.setFrontID(mainID);
//                            fusionPattern.setBackID(subID);
//                            fusionPattern.setOverlapBase(overlapBase);
//                        }else if(mainEndIndex>subEndIndex){
//                            // main is back part
//                            long overlapBase = (subIniIndex - subEndIndex);
//                            frontCode = (subChrID<<28)+subEndMatchPos;
//                            backCode = (mainChrID<<28)+mainIniMatchPos;
//                            fusionPattern.setFrontID(subID);
//                            fusionPattern.setBackID(mainID);
//                            fusionPattern.setOverlapBase(overlapBase);
//                        }else if(subEndIndex==mainEndIndex){
//                            continue;
//                        } 
                        continue;
                    }else{
                        continue;
                    }
                    
                    /**
                     * Group fusion pattern
                     * Store in Map
                     */

                    if(fusionGroupMap.containsKey(frontCode)){
                        Map<Long,ArrayList<LASTFusionPattern>> layerTwo = fusionGroupMap.get(frontCode);
//                        if(frontCode == 1648306277 && backCode == 1192524689){
//                            System.out.println();
//                        }
                        
                        if(layerTwo.containsKey(backCode)){
                            ArrayList<LASTFusionPattern> listFusionPattern = layerTwo.get(backCode);
                            listFusionPattern.add(fusionPattern);
                            layerTwo.put(backCode, listFusionPattern);
                        }else{
                            ArrayList<LASTFusionPattern> listFusionPattern = new ArrayList();
                            listFusionPattern.add(fusionPattern);
                            layerTwo.put(backCode, listFusionPattern);
                        }
                        
                        fusionGroupMap.put(frontCode, layerTwo);
                    }else{
                        Map<Long,ArrayList<LASTFusionPattern>> layerTwo = new LinkedHashMap();
                        ArrayList<LASTFusionPattern> listFusionPattern = new ArrayList();
                        listFusionPattern.add(fusionPattern);
                        layerTwo.put(backCode, listFusionPattern);
                        fusionGroupMap.put(frontCode, layerTwo);
                    }
                }
            }   
        }
        
        
        /**
         * Write FusionGroup to File
         */
//        String[] dummyfilename = intersectFile.split("/");       
//        String sfilename = dummyfilename[dummyfilename.length-1].split("\\.")[0];
        String saveFileName = intersectFile.split("\\.")[0] + "_Fusion.txt";
        FileWriter writer;        
        /**
         * Check File existing
         */
        ArrayList<Integer> countRead = new ArrayList();
        File f = new File(saveFileName); //File object        
        if(f.exists()){
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(saveFileName,true);
        }else{
//            ps = new PrintStream(filename);
            writer = new FileWriter(saveFileName);
        }
        
        int numGroupCount = 0;

        for(Map.Entry<Long,Map<Long,ArrayList<LASTFusionPattern>>> entry : fusionGroupMap.entrySet()){
            long fCode = entry.getKey();
            Map<Long,ArrayList<LASTFusionPattern>> lTwo = entry.getValue();
            
            for(Map.Entry<Long,ArrayList<LASTFusionPattern>> entry2 : lTwo.entrySet()){
                long bCode = entry2.getKey();
                ArrayList<LASTFusionPattern> patternList = entry2.getValue();
                int coverage = patternList.size();
                
                if(coverage >= threshold){
                    writer.write("Group "+(++numGroupCount)+"\tCoverage "+coverage+"\n");
                    for(int i=0;i<patternList.size();i++){
                        LASTFusionPattern fusionPattern = patternList.get(i);
                        int fID = fusionPattern.getFrontID();
                        int bID = fusionPattern.getBackID();
                        long overlapBase = fusionPattern.getOverlapBase();
                        
                        LASTResult fLAST = listLASTResult.get(fID);
                        LASTResult bLAST = listLASTResult.get(bID);
                        
                        writer.write(fLAST.getScore()+"\t"+fLAST.getRefChr()+"\t"+fLAST.getRefPos()+"\t"+fLAST.getRefNumBase()+"\t"+fLAST.getRefStrand()+"\t"+fLAST.getRefSize()+"\t"+fLAST.getSampleName()+"\t"+fLAST.getSamplePos()+"\t"+fLAST.getSampleNumBase()+"\t"+fLAST.getSampleStrand()+"\t"+fLAST.getSampleSize()+"\t"+fLAST.getSampleSequenceMatch()
                                +"\t"+bLAST.getScore()+"\t"+bLAST.getRefChr()+"\t"+bLAST.getRefPos()+"\t"+bLAST.getRefNumBase()+"\t"+bLAST.getRefStrand()+"\t"+bLAST.getRefSize()+"\t"+bLAST.getSampleName()+"\t"+bLAST.getSamplePos()+"\t"+bLAST.getSampleNumBase()+"\t"+bLAST.getSampleStrand()+"\t"+bLAST.getSampleSize()+"\t"+bLAST.getSampleSequenceMatch()
                                +"\t"+overlapBase+"\n");   
                    }
                    writer.write("\n");
                }
            }
        }
        
        writer.flush();
        writer.close();        
    }
    
    public static ArrayList<LASTResult> readLASTResult(String fileName) throws IOException{
        Charset charset = Charset.forName("US-ASCII");
        
        ArrayList<LASTResult> listLASTResult = new ArrayList();
        
        int score=0;
        String EG2="";
        String E="";
        
        String refChr="";
        long refPos=0;
        int refNumBase=0;
        String refStrand="";
        int refSize=0;
        String refSequenceMatch="";
        long refSequenceCode=0;
        
        String sampleName="";
        long samplePos=0;
        int sampleNumBase=0;
        String sampleStrand="";
        int sampleSize=0;
        String sampleSequenceMatch="";
        long sampleSequenceCode=0;
        
        //String[] ddSS = filename.split(".");
//        String saveFileName = inLASTtResultFile.split("\\.")[0] + "_LASTCountSample.txt";
        Path path = Paths.get(fileName);
        Map<String,Integer> mapCountSample = new LinkedHashMap();
        int dummy = 1;
        String data = "";
        ArrayList<String> inData = new ArrayList();    
        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            
            while ((data = reader.readLine()) != null) {
                if(data.isEmpty() == false){ 
                    if(data.charAt(0) == 's' && dummy == 1){
                        dummy++;
                        String[] splitData = data.split("\\s+");
                        refChr = splitData[1];
                        refPos = Long.parseLong(splitData[2]);
                        refNumBase = Integer.parseInt(splitData[3]);
                        refStrand = splitData[4];
                        refSize = Integer.parseInt(splitData[5]);
                        refSequenceMatch = splitData[6];
                        refSequenceCode = SequenceUtil.encodeMer(refSequenceMatch, 18);
                        
                    }else if(data.charAt(0) == 's' && dummy == 2) {
                        dummy = 1;
                        
                        String[] splitData = data.split("\\s+");
                        sampleName = splitData[1];                        
                        samplePos = Long.parseLong(splitData[2]);
                        sampleNumBase = Integer.parseInt(splitData[3]);
                        sampleStrand = splitData[4];
                        sampleSize = Integer.parseInt(splitData[5]);
                        sampleSequenceMatch = splitData[6];
                        sampleSequenceCode = SequenceUtil.encodeMer(sampleSequenceMatch, 18); 

                        if(mapCountSample.containsKey(sampleName)){
                            int amount = mapCountSample.get(sampleName);
                            mapCountSample.put(sampleName, ++amount);
                        }else{
                            mapCountSample.put(sampleName, 1);
                        }
                        
                        LASTResult lastRes = new LASTResult();
                        lastRes.addRefData(refChr,refPos,refNumBase,refStrand,refSize,refSequenceMatch,refSequenceCode);
                        lastRes.addSampleData(sampleName,samplePos,sampleNumBase,sampleStrand,sampleSize,sampleSequenceMatch,sampleSequenceCode);
                        lastRes.addScoreData(score,EG2,E);
                        listLASTResult.add(lastRes);
                        
                    }else if(data.charAt(0)=='a'){
                        String[] splitData = data.split("\\s+");
                        score = Integer.parseInt(splitData[1].split("=")[1]);
                        EG2 = splitData[2].split("=")[1];
                        E = splitData[3].split("=")[1];
                    }                
                }
            }                
        }
        
        return listLASTResult;
    }
}
