/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ClusterGroup;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import static java.lang.Long.min;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author worawich
 */
public class Clustering {
    
    public static ArrayList<ClusterGroup> clusteringGroup(AlignmentResultRead inAlnRead, double threshold){
        ArrayList<ClusterGroup> listGroup = new ArrayList();
        ArrayList checkList = new ArrayList();
        inAlnRead.createGroupCharacteristic(threshold); // create significant data for clustering purpose
        ArrayList<ShortgunSequence> listSS = inAlnRead.getResult();
        ClusterGroup group = new ClusterGroup();
        for(int i=0;i<listSS.size();i++){
            ShortgunSequence mainDummySS = listSS.get(i);
            ArrayList mainOutGroup = mainDummySS.getOutGroup();
            
            if(checkList.contains(i)){
                
            }else{
                group = new ClusterGroup();
                group.addShortgunRead(mainDummySS);
                checkList.add(i);
                
                for(int j=0;j<listSS.size();j++){
                    
                    if(j!=i){
                        ShortgunSequence subDummySS = listSS.get(j);
                        ArrayList subOutGroup = subDummySS.getOutGroup();
                        
                        if(mainOutGroup.equals(subOutGroup)){
                            
                            
                            int strandCheck = checkStrand(mainDummySS,subDummySS);
                            
                            if(strandCheck == 1){
                                group.addShortgunRead(subDummySS);
                                checkList.add(j);
                            }
                            
                        }
                    }
                }
                listGroup.add(group);
            }     
        }
        return listGroup;   
    }
    
    public static Map<Long,ArrayList<String>> createChrMatchMap(AlignmentResultRead inAlnRead){
        /* 
            Create chromosome table (P' nung guide)
        */
        
        ArrayList<ShortgunSequence> listSS = inAlnRead.getResult();
        Map<Long,ArrayList<String>> classMap = new HashMap(); 
        
        for(int i=0;i<listSS.size();i++){
            ShortgunSequence dummySS = listSS.get(i);
            
            ArrayList chrMatch = dummySS.getListChrMatch();
            
            if(i==0){
                for(long j=1;j<25;j++){
                    if(chrMatch.contains(j)){
                        ArrayList<String> dummyReadList = new ArrayList();
                        dummyReadList.add(dummySS.getReadName());
                        classMap.put(j, dummyReadList);  
                    }else{
                        ArrayList<String> dummyReadList = new ArrayList();
                        classMap.put(j, dummyReadList); 
                    }
                }
                
            }else{
                
                for(long j=1;j<25;j++){
                    if(chrMatch.contains(j)){
                        ArrayList<String> dummyReadList = classMap.get(j);
                        dummyReadList.add(dummySS.getReadName());
                        classMap.put(j, dummyReadList);

                    }

                }
            }
            
        }
        
        return classMap;
    }
    
    public static Map<Integer,ArrayList<String>> filterClusterGroupLocalAlignment(ArrayList<Map<Integer,ArrayList<String>>> indata, int minCoverage){
        Map dummyMap = new HashMap();
        Map<Integer,ArrayList<String>> mainGroupMap = new HashMap();
        int numGroup = 0;
        
        for(int i=0 ; i<indata.size() ;i++){
            dummyMap = indata.get(i);
            
            if(dummyMap.isEmpty()!=true){
                
                Set dummySet = dummyMap.keySet();
                Iterator dummyIter = dummySet.iterator();
                while(dummyIter.hasNext()){
                    Integer dummyKey = (Integer) dummyIter.next();
                    ArrayList<String> group = (ArrayList<String>) dummyMap.get(dummyKey);
                    
                    if(group.size()>minCoverage && mainGroupMap.containsValue(group)!=true){
                        numGroup++;
                        mainGroupMap.put(numGroup, group);
                    }
                }
            }
        }
        
        
        
        return mainGroupMap;
    }
    
    public static void writeLocalAlignmentInFile(Map<Integer,ArrayList<String>> inData, String path, String filename) throws FileNotFoundException{
        
        PrintStream ps = new PrintStream(path+filename+".txt");            // Create file object
        
        for(Map.Entry<Integer,ArrayList<String>> entry : inData.entrySet()){
            Integer numGroup = entry.getKey();
            ArrayList<String> group = entry.getValue();
            
            ps.print("Group" + numGroup + ":");
            
            for(int i=0;i<group.size();i++){
                ps.print("\t"+group.get(i));
            }
            ps.println();
        }
    }
    
    public static ArrayList<ClusterGroup> clusteringGroupV2(AlignmentResultRead inAlnRead, double threshold){
        /* New implementation cluster without create ClusterGroup but use Map instead */
        /*  In process  */
        ArrayList<ClusterGroup> listGroup = new ArrayList();
        ArrayList checkList = new ArrayList();
        inAlnRead.createGroupCharacteristic(threshold); // create significant data for clustering purpose
        ArrayList<ShortgunSequence> listSS = inAlnRead.getResult();
        ClusterGroup group = new ClusterGroup();
        for(int i=0;i<listSS.size();i++){
            ShortgunSequence mainDummySS = listSS.get(i);
            ArrayList mainOutGroup = mainDummySS.getOutGroup();
            
            if(checkList.contains(i)){
                
            }else{
                group = new ClusterGroup();
                group.addShortgunRead(mainDummySS);
                checkList.add(i);
                
                for(int j=0;j<listSS.size();j++){
                    
                    if(j!=i){
                        ShortgunSequence subDummySS = listSS.get(j);
                        ArrayList subOutGroup = subDummySS.getOutGroup();
                        
                        if(mainOutGroup.equals(subOutGroup)){
                            
                            
                            int strandCheck = checkStrand(mainDummySS,subDummySS);
                            
                            if(strandCheck == 1){
                                group.addShortgunRead(subDummySS);
                                checkList.add(j);
                            }
                            
                        }
                    }
                }
                listGroup.add(group);
            }     
        }
        return listGroup;   
    }
    
    public static int checkStrand(ShortgunSequence mainSS, ShortgunSequence subSS){
        // Check strand of match result protect from miss group classified in case that both read are align in same chromosome and align position but different strand
        // similarity = 0 mean not the same || similarity = 1 mean it the same
        int check = 0;
        int similarity = 0;
        int sizeMatchMainSS = mainSS.getListChrMatch().size();
        int sizeMatchSubSS = subSS.getListChrMatch().size();
        
        ArrayList matchMainChr = mainSS.getListChrMatch();
        ArrayList matchMainStrand = mainSS.getListStrand();
        ArrayList matchSubChr = subSS.getListChrMatch();
        ArrayList matchSubStrand = subSS.getListStrand();
        
        if(sizeMatchMainSS == sizeMatchSubSS){
            for(int i=0;i<sizeMatchMainSS;i++){
                int mainChr = (int)matchMainChr.get(i);
                String mainStrand = (String)matchMainStrand.get(i);
                    
                for(int j=0;j<sizeMatchMainSS;j++){
                    int subChr = (int)matchSubChr.get(j);
                    String subStrand = (String)matchSubStrand.get(j);

                    if(mainChr == subChr & mainStrand.equals(subStrand)){
                        check = check+1;
                    }   
                }     
            }
        }
        
        
        if(check == sizeMatchMainSS){
            similarity = 1;
        }
        
        return similarity;
    }
    
    public static void createColorArray(AlignmentResultRead inRes, int rLen, int merLen){
        /**
         *  Use for create colorArray 
         *  please make sure that your inRes must implement iniIndex in its
         *  
         *  Suitable for version 3 data structure (data structure that has iniIdx in its)
         */
        
        
        int iniIdx;
        int lastIdx;
        int numMatch;
        int numMer = (int) (rLen-merLen)+1;
        int[] matchArray = new int[numMer];
        ArrayList<ShortgunSequence> listRead = inRes.getResult();
        ArrayList<int[]> listMatchArray = new ArrayList();
        Map<String,String[]> colorPatternMap = new HashMap();
        
        for(int i =0;i<listRead.size();i++){ // Loop all Read
            ShortgunSequence dummySS = listRead.get(i);            
            
            ArrayList<Integer> listChr = dummySS.getListChrMatch();                // all of Array below has same size (size indicate the number of match pattern)
            ArrayList<Long> listPos = dummySS.getListPosMatch();
            ArrayList<String> listStrand = dummySS.getListStrand();
            ArrayList<Integer> listNumMatch = dummySS.getListNumMatch();
            ArrayList<Integer> listIniIdx = dummySS.getListIniIdx();
            long numPattern = listChr.size();
            int[] iniIdxArray = new int[(int)numPattern];
            int[] lastIdxArray = new int[(int)numPattern];
            /**
             * Start create match Array of each pattern of this read
             *      integer array size numMer
             *  code 1 = match ; code 0 = not match
             */
            listMatchArray = new ArrayList();
            for(int numP=0;numP<numPattern;numP++){             // loop over match pattern
                matchArray = new int[numMer];
               
                numMatch = listNumMatch.get(numP);
                String strand = listStrand.get(numP);
                
                if(strand.equals("-")){                          // check for strand type (the iniIndx of type - must have recalculate with the equation in this case
                    int dummyIniIdx = listIniIdx.get(numP);
                    iniIdx = rLen - (dummyIniIdx + ((merLen+numMatch)-1));
                    lastIdx = (iniIdx + numMatch)-1;
//                    System.out.println("iniIdx check" + iniIdx);
                }else{
                    iniIdx = listIniIdx.get(numP);
                    lastIdx = (iniIdx + numMatch)-1;
                }
                /**
                 * add initial index and last index of each pattern for future use
                 */
                iniIdxArray[numP] = iniIdx;
                lastIdxArray[numP] = lastIdx;
                
                
                for(int merIdx = (int)iniIdx; merIdx<(iniIdx+numMatch); merIdx++){
                    matchArray[merIdx] = 1;      
                }
                listMatchArray.add(matchArray);                 // At the end, this array should be the same size as number of match pattern         
            }
            
            /**
             * Start create color array of this read
             * 
             */
            byte[] colorArray = new byte[numMer];                   // We reconsider to store in byte(8bit) not char(16bit) to save memory cousumpsion
            for(int index=0;index<numMer;index++){                  // Loop over index of possible mer (100-18)+1 = 83 (max index)
                               
                ArrayList<Integer> chrCheckList = new ArrayList();
                for(int p=0;p<listMatchArray.size();p++){           // Loop over listMatchArray size it is seemlessly as we loop over match pattern
                    int[] dummyMatchArray = listMatchArray.get(p);
                    
                    if(dummyMatchArray[index]==1){
                        chrCheckList.add(listChr.get(p));
                    }  
                }
                
                if(chrCheckList.isEmpty()){
                    colorArray[index] = 0;    // chrCheckList is empty, this mean no match at this index. we assign 0 to colorArray
                }else{
                    /**
                     *  check duplicate from chrCheckList
                     *  and assign color string to each mer index
                     */
                    Set<Integer> chrCheckSet = new HashSet<Integer>(chrCheckList);
                    if(chrCheckSet.size()<chrCheckList.size()){         // check case if this true that mean chrCheckList has duplicate element in it
                        if(chrCheckSet.size() == 1){    // Has duplicate and size is 1, this mean it match only one chr but various position
                            colorArray[index] = 3;
                        }else{                          // Has duplicate and size is more than 1, this mean it match at same chr and other chr
                            colorArray[index] = 4;
                        }
                    }else{
                        if(chrCheckList.size()>1){      // Has no duplicate and size is more than 1, this mean it match at different chr
                            colorArray[index] = 2;
                        }else if(chrCheckList.size()==1){   // Has no duplicate and size is 1, this mean it unique
                            colorArray[index] = 1;
                        }
                    }
                }
 
            }
           
            /**
             *  Store colorArray in Map which has key=read name and value=colorArray (has size equal to number of possible mer)
             */
            //colorPatternMap.put(dummySS.getReadName(), colorArray);
            dummySS.addReadLength(rLen);
            dummySS.addMerLength(merLen);
            dummySS.addColorArray(colorArray);      // store colorArray in Shortgun sequence (not sure with this one hope that store in the AlignmentResultRead that we feed as input when this function finish)
            
        }
        
    }
    
    public static ArrayList<ClusterGroup> clusterFromFile(String filename, int maxBaseDiff, int minCoverage) throws IOException{
        /**
         * Suitable for linux sort format only
         */
        
        ClusterGroup group = new ClusterGroup();
        ArrayList<ClusterGroup> listGroup = new ArrayList();
        Charset charset = Charset.forName("US-ASCII");
        //String[] ddSS = filename.split(".");
        String saveFileName = filename.split("\\.")[0] + "_ClusterGroup.txt";
        Path path = Paths.get(filename);

        StringBuffer seq = new StringBuffer();
        ArrayList<String> inData = new ArrayList();    
        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;
            String[] data = null;
            long oldIniPos = 0;
            long diff = 0;
            byte oldNumChr = 0;
            int count = 0;
            
            System.out.println("reading");
            while ((line = reader.readLine()) != null) {
                
                inData.add(line);
                count++;
                if(count%1000000==0){
                    System.out.println(count + " line past");
                    //System.out.println("Recent chromosome: " + numChr);
                }       
                
            }
//            writeClusterGroupToFile(filename,listGroup);
        }
        
        System.out.println("Done reading");
        System.out.println("Grouping");
        
        long oldIniPos = 0;
        long diff = 0;
        byte oldNumChr = 0;
        int count = 0;
        String data = null;
        String[] splitData = new String[10]; 
        for(int i=0;i<inData.size();i++){                       // loop get each data
            data = inData.get(i);
            byte numChr = 0;
            long iniPos = 0;
            String readName = null;
            int field = 0;
                
            int lastIndex = 0;
            for(int index=0;index<data.length();index++){       // loop split data 
                char c = data.charAt(index);

                if(c == ','){
                    field++;
                    switch(field)
                    {
                        case 1:
                            numChr = Byte.parseByte(data.substring(lastIndex, index));
                            break;
                        case 2:
                            iniPos = Long.parseLong(data.substring(lastIndex, index));
                            break;
                        case 9:
                            readName = data.substring(index+1,data.length());         // get the last string                            
                            break;
                    }                    

                    lastIndex = index + 1;  
                }
            }
                
            diff = iniPos - oldIniPos;

            if(oldNumChr == numChr){            // same Group check first criteria
                if(diff <= maxBaseDiff){        // same Group check second criteria
                    if(group.getListReadname().contains(readName)!=true){
                        group.addReadName(readName);
                        group.addChromosomeNumber(numChr);
                        group.addIniPos(iniPos);
                    }

                }else{

                    if(group.getNumMember() > minCoverage){
                        listGroup.add(group);
                    }else{
                        group = null;
                        System.gc();
                    }

                    group = new ClusterGroup();

                    group.addReadName(readName);
                    group.addChromosomeNumber(numChr);
                    group.addIniPos(iniPos);                 
                }
            }else{

                if(group.getNumMember() > minCoverage){
                    listGroup.add(group);

                }else{
                    group = null;
                    System.gc();
                }

                group = new ClusterGroup();

                group.addReadName(readName);
                group.addChromosomeNumber(numChr);
                group.addIniPos(iniPos);
            }

            oldNumChr = numChr;
            oldIniPos = iniPos;
            count++;
            if(count%1000000==0){
                System.out.println(count + " match pattern past");
                System.out.println("Recent chromosome: " + numChr);
            }
                
            
                
//                byte numChr = Byte.parseByte(data[0]);
//                long iniPos = Long.parseLong(data[1]);
////                long lastPos = Long.parseLong(data[2]);
////                byte numG = Byte.parseByte(data[3]);
////                byte numY = Byte.parseByte(data[4]);
////                byte numO = Byte.parseByte(data[5]);
////                byte numR = Byte.parseByte(data[6]);
////                String strand = data[7];
////                byte iniIdx = Byte.parseByte(data[8]);
//                String readName = data[9];
//                diff = iniPos - oldIniPos;
//                
//                if(oldNumChr == numChr){            // same Group check first criteria
//                    if(diff <= maxBaseDiff){        // same Group check second criteria
//                        if(group.getListReadname().contains(readName)!=true){
//                            group.addReadName(readName);
//                            group.addChromosomeNumber(numChr);
//                            group.addIniPos(iniPos);
////                            group.addLastPos(lastPos);
////                            group.addNumGreen(numG);
////                            group.addNumYellow(numY);
////                            group.addNumOrange(numO);
////                            group.addNumRed(numR);
////                            group.addStrand(strand);
////                            group.addIniIndex(iniIdx);
//                        }
//                        
//                    }else{
//                        
//                        if(group.getNumMember() > minCoverage){
//                            listGroup.add(group);
//                        }else{
//                            group = null;
//                            System.gc();
//                        }
//                            
//                        group = new ClusterGroup();
//                        
//                        group.addReadName(readName);
//                        group.addChromosomeNumber(numChr);
//                        group.addIniPos(iniPos);
////                        group.addLastPos(lastPos);
////                        group.addNumGreen(numG);
////                        group.addNumYellow(numY);
////                        group.addNumOrange(numO);
////                        group.addNumRed(numR);
////                        group.addStrand(strand);
////                        group.addIniIndex(iniIdx);                       
//                    }
//                }else{
//                    
//                    if(group.getNumMember() > minCoverage){
//                        listGroup.add(group);
//                        
//                    }else{
//                        group = null;
//                        System.gc();
//                    }
//                    
//                    group = new ClusterGroup();
//
//                    group.addReadName(readName);
//                    group.addChromosomeNumber(numChr);
//                    group.addIniPos(iniPos);
////                    group.addLastPos(lastPos);
////                    group.addNumGreen(numG);
////                    group.addNumYellow(numY);
////                    group.addNumOrange(numO);
////                    group.addNumRed(numR);
////                    group.addStrand(strand);
////                    group.addIniIndex(iniIdx);
//                }
//                
//                oldNumChr = numChr;
//                oldIniPos = iniPos;
//                count++;
//                if(count%1000000==0){
//                    System.out.println(count + " match pattern past");
//                    System.out.println("Recent chromosome: " + numChr);
//                }  
        }
        
        System.out.println("Done Grouping");
        System.out.println("Writing");
        writeClusterGroupToFile(filename,listGroup);
        System.out.println("Done Writing");
        return listGroup;
    }
    
    public static ArrayList<ClusterGroup> clusterFromFileBinarySearch(String filename, int maxBaseDiff, int minCoverage, int significantCluster) throws IOException{
    /**
     * Suitable for linux sort format only
     */

    ClusterGroup group = new ClusterGroup();
    ArrayList<ClusterGroup> listGroup = new ArrayList();
    Charset charset = Charset.forName("US-ASCII");
    //String[] ddSS = filename.split(".");
    String saveFileName = filename.split("\\.")[0] + "_ClusterGroup.txt";
    Path path = Paths.get(filename);

    StringBuffer seq = new StringBuffer();
    ArrayList<String> inData = new ArrayList();    
    try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
        String line = null;
        String[] data = null;
        long oldIniPos = 0;
        long diff = 0;
        byte oldNumChr = 0;
        int count = 0;

        System.out.println("reading");
        while ((line = reader.readLine()) != null) {

            inData.add(line);
            count++;
            if(count%1000000==0){
                System.out.println(count + " line past");
                //System.out.println("Recent chromosome: " + numChr);
            }       

        }
//            writeClusterGroupToFile(filename,listGroup);
    }

    System.out.println("Done reading");
    System.out.println("Grouping");

    listGroup = clusterBinaryImplement(inData,maxBaseDiff,minCoverage,significantCluster);
    
    inData.clear();
    System.gc();
    
    System.out.println("Done Grouping");
    System.out.println("Writing");
    writeClusterGroupToFile(filename.split("\\.")[0] + "_ClusterGroup_first_Phase.txt",listGroup);
    System.out.println("Done Writing");
    System.out.println("Start second Grouping");
    System.out.println("extract and save signifcant group");
    ArrayList<ClusterGroup> sigGroup = clusterPhaseTwo_groupBySignificant(listGroup);
    writeClusterGroupToFile(filename.split("\\.")[0] + "_ClusterGroup_second_Phase_sig.txt",sigGroup);
    //ArrayList<ClusterGroup> secondListGroup = createGroupIntersect(listGroup);
    System.out.println("Done extract and save significant Group");
    System.out.println("Start extract common peak");
    Map<String,ArrayList<Integer>> readMap = clusterPhaseTwo_createReadMap(sigGroup);
    ArrayList<ClusterGroup> lastListGroup = clusterPhaseTwo_extractCommon(sigGroup,readMap);
    System.out.println("Done extract common peak");
    System.out.println("Save extract common peak");
    writeClusterGroupToFile(filename.split("\\.")[0] + "_ClusterGroup_second_Phase_CommonPeaks.txt",lastListGroup);
    System.out.println("Done save extract common peak");
//    System.out.println("Start group intersect");
//    ArrayList<ClusterGroup> secondListGroup = createGroupIntersect(sigGroup);
//    System.out.println("Done group intersect");
//    System.out.println("Writing");
//    writeClusterGroupToFile(saveFileName,secondListGroup);
//    System.out.println("Done Writing");
    return lastListGroup;
}

    
    public static void writeClusterGroupToFile(String filename,ArrayList<ClusterGroup> input) throws IOException{
        ArrayList<String> readNameList;                             // use for local alignment and other stuff
        ArrayList<Byte> listChr;
        ArrayList<Long> listIniPos;
        ArrayList<Long> listLastPos;
        ArrayList<Byte> listNumG;
        ArrayList<Byte> listNumY;
        ArrayList<Byte> listNumO;
        ArrayList<Byte> listNumR;       
        ArrayList<String> listStrand;       
        ArrayList<Byte> listIniIndex;
        
        
        File file = new File(filename);
        
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
        
        for(int i=0; i<input.size();i++){                   // loop over each group
            ClusterGroup dummyGroup = input.get(i);
            
            readNameList = dummyGroup.getListReadname();
            listChr = dummyGroup.getListChromosome();
            listIniPos = dummyGroup.getListIniPos();
            listLastPos = dummyGroup.getListLastPos();
            listNumG = dummyGroup.getListNumGreen();
            listNumY = dummyGroup.getListNumYellow();
            listNumO = dummyGroup.getListNumOrange();
            listNumR = dummyGroup.getListNumRed();
            listStrand = dummyGroup.getListStrand();
            listIniIndex = dummyGroup.getListIniIndex();            
            
            writer.write(">Group "+i);
            writer.write("\n");
            for(int j=0; j<dummyGroup.getNumMember();j++){
                if(listChr.isEmpty()){
                    writer.write(String.format("%s",readNameList.get(j)));
                    writer.write("\n");
                }else{
                    writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s", listChr.get(j),listIniPos.get(j),listLastPos.get(j),listNumG.get(j),listNumY.get(j),listNumO.get(j),listNumR.get(j),listStrand.get(j),listIniIndex.get(j),readNameList.get(j)));
                    writer.write("\n");
                }   
            }
            
        }
        
        writer.flush();
        writer.close();
    }
    
    public static ArrayList<ClusterGroup> clusterBinaryImplement(ArrayList<String> inData,int maxBaseDiff, int minCoverage, int threshold){     // thrshold use in criteria for extract significant read from bunct of read
       
        ArrayList<String> listData = inData; 
        int dataIdx = 0;
        int lowIdx = 0;
        int highIdx = listData.size()-1;
        int lastIdx = 0;
        ClusterGroup group = new ClusterGroup();
        ArrayList<ClusterGroup> listGroup = new ArrayList();
        int count = 0;
        while(true){

            String compareData = listData.get(lowIdx);
            /***    Extract data    ****/
            String[] data = compareData.split(",");
            byte numChr = Byte.parseByte(data[0]);
            long iniPos = Long.parseLong(data[1]);
            long lastPos = Long.parseLong(data[2]);
            byte numG = Byte.parseByte(data[3]);
            byte numY = Byte.parseByte(data[4]);
            byte numO = Byte.parseByte(data[5]);
            byte numR = Byte.parseByte(data[6]);
            String strand = data[7];
            byte iniIdx = Byte.parseByte(data[8]);
            String readName = data[9];
            /*******/
            
            /***********    Do binary search ***********/
            int dummyHighIdx = binarySearch(lowIdx, highIdx, listData, numChr, iniPos, maxBaseDiff);
            
            /*********  Do scanning down 200 line **********/
            
            if(dummyHighIdx > 0){
                for(int i=dummyHighIdx+1;i<=min(highIdx,dummyHighIdx+200);i++){
                     
                    String[] refData = listData.get(i).split(",");
                    byte refNumChr = Byte.parseByte(refData[0]);
                    long refIniPos = Long.parseLong(refData[1]);
                    long refLastPos = Long.parseLong(refData[2]);
                    byte refNumG = Byte.parseByte(refData[3]);
                    byte refNumY = Byte.parseByte(refData[4]);
                    byte refNumO = Byte.parseByte(refData[5]);
                    byte refNumR = Byte.parseByte(refData[6]);
                    String RefStrand = refData[7];
                    byte refIniIdx = Byte.parseByte(refData[8]);
                    String refReadName = refData[9];
                    
                    long diff = iniPos - refIniPos;
                    
                    if(diff>maxBaseDiff||diff<maxBaseDiff*-1){
                        lastIdx = i-1;
                        break;
                    }else{
                        lastIdx = i;
                    }                         
                }
                
                /******* create cluster group and add group *******/
                 
                group = new ClusterGroup();
                for(int idx=lowIdx;idx<=lastIdx;idx++){
                    String[] realData = listData.get(idx).split(",");
                    numChr = Byte.parseByte(realData[0]);
                    iniPos = Long.parseLong(realData[1]);
                    lastPos = Long.parseLong(realData[2]);
                    numG = Byte.parseByte(realData[3]);
                    numY = Byte.parseByte(realData[4]);
                    numO = Byte.parseByte(realData[5]);
                    numR = Byte.parseByte(realData[6]);
                    strand = realData[7];
                    iniIdx = Byte.parseByte(realData[8]);
                    readName = realData[9];
                    if(group.getListReadname().contains(readName)!=true){
                        group.addReadName(readName);
                        group.addChromosomeNumber(numChr);
                        group.addIniPos(iniPos);

                        /* Additional information */
                        group.addIniIndex(iniIdx);
                        group.addLastPos(lastPos);
                        group.addNumGreen(numG);
                        group.addNumOrange(numO);
                        group.addNumRed(numR);
                        group.addNumYellow(numY);
                        group.addStrand(strand);



                        int numMerMatch = numG+numY;

                        if(numMerMatch >= threshold){
                            group.significantSetTrue();
                            group.addHighlightRead(readName);
                        }
                    } 
                }
                
                /***********/
            }else{
                // no coverage ignore create new group
                lastIdx = dummyHighIdx;
//                group = new ClusterGroup();
//                group.addReadName(readName);
//                group.addChromosomeNumber(numChr);
//                group.addIniPos(iniPos);  
            }
            
            if(group.getNumMember() > minCoverage){       // check case for add only high coverage
                listGroup.add(group);
            }
            
            lowIdx = lastIdx+1;
            if(lastIdx==highIdx||lowIdx==highIdx){
                break;
            }
            count++;
            if(count%100000==0){
                System.out.println(count + " round past");
                System.out.println("Current next index is : "+lowIdx);
                //System.out.println("Recent chromosome: " + numChr);
            }  
        }
        
        listData.clear();
        System.gc();
        
        return listGroup;
    }
    
    public static int binarySearch(int lowIdx, int highIdx, ArrayList<String> listRefData,byte inChr,long inIniPos,int threshold){
        
        while(lowIdx <= highIdx){
            
            int middle = (lowIdx + highIdx)/2;
            if(middle == lowIdx){
                //middle++;
            }
            /***    Extract data    ****/
            String[] refData = listRefData.get(middle).split(",");       
            byte numChr = Byte.parseByte(refData[0]);
            long iniPos = Long.parseLong(refData[1]);
            long lastPos = Long.parseLong(refData[2]);
            byte numG = Byte.parseByte(refData[3]);
            byte numY = Byte.parseByte(refData[4]);
            byte numO = Byte.parseByte(refData[5]);
            byte numR = Byte.parseByte(refData[6]);
            String strand = refData[7];
            byte iniIdx = Byte.parseByte(refData[8]);
            String readName = refData[9];
            /*******/
            
            long diff = inIniPos - iniPos;
            
            
            if(numChr == inChr){
                if(diff>=(threshold*-1)&&diff<=threshold){
                    return middle;
                }else if(diff > threshold){
                    lowIdx = middle;
                    middle = (lowIdx + highIdx)/2;
                }else if(diff < (threshold*-1)){
                    highIdx = middle;
                    middle = (lowIdx + highIdx)/2;
                }

            }else if(numChr > inChr){
                highIdx = middle;
                middle = (lowIdx + highIdx)/2;
            }else if(numChr < inChr){
                lowIdx = middle;
                middle = (lowIdx + highIdx)/2;            
            }
        }
        
        return -1;
    }
    
    public static ArrayList<ClusterGroup> createGroupIntersect(ArrayList<ClusterGroup> inListGroup){

        /**
         *  First part (create characteristic)
         */
        
        Map<String,boolean[]> groupMap = new HashMap();
        boolean[] groupChar = new boolean[inListGroup.size()];                  // initialize boolean array as group characteristic for each read
        ArrayList<String> checkList = new ArrayList();
        int count = 0;
        
        for(int mainG=0;mainG<inListGroup.size();mainG++){                      // Loop over each group
            ClusterGroup dummyGroup = inListGroup.get(mainG);
            ArrayList<String> dummyListRead = dummyGroup.getListReadname();
            
            for(int numMainRead=0;numMainRead<dummyListRead.size();numMainRead++){  // Loop over each read in selected group
                String mainRead = dummyListRead.get(numMainRead);
                
                if(checkList.contains(mainRead)!=true){
                    groupChar = new boolean[inListGroup.size()];
                    for(int subG=0;subG<inListGroup.size();subG++){                     // Loop over all group to check containment of selected read
                        ClusterGroup subGroup = inListGroup.get(subG);
                        ArrayList<String> subListRead = subGroup.getListReadname();

                        if(subListRead.contains(mainRead)){
                            groupChar[subG]=true;
                        }else{
                            groupChar[subG]=false;
                        }
                    }
                    groupMap.put(mainRead, groupChar);
                    checkList.add(mainRead);
                }   
            }
            
            count++;
              if(count%50000==0){
                System.out.println("Create characteristic "+count + " round past");
            }
        }
        
        
        /**
         * Second part (reGroup with intersect characteristic)
         */
        ArrayList<ClusterGroup> listGroup = new ArrayList();
        checkList = new ArrayList();
        count = 0;
        for(Map.Entry<String,boolean[]> mainEntry : groupMap.entrySet()){       // main Loop for each read
            String mainKey = mainEntry.getKey();
            boolean[] mainValue = mainEntry.getValue();
            ClusterGroup group = new ClusterGroup();
            
            if(checkList.contains(mainKey)!=true){
                group.addReadName(mainKey);
                checkList.add(mainKey);
                
                for(Map.Entry<String,boolean[]> subEntry : groupMap.entrySet()){    // minor loop for check intersect
                    String subKey = subEntry.getKey();
                    boolean[] subValue = subEntry.getValue();
                    if(checkList.contains(subKey)!=true){
                        if(Arrays.equals(mainValue, subValue)){
                            group.addReadName(subKey);
                        }
                    }   
                }
            }
            
            listGroup.add(group);
            count++;
            if(count%50000==0){
                System.out.println("Regroup with intersect characteristic "+count + " round past");
            }
        }         
        return listGroup;
    }
    
    public static void clusterPhaseTwo(ArrayList<ClusterGroup> inListGroup){
        
        Map<String,ArrayList<Integer>> readNameMap = clusterPhaseTwo_createReadMap(inListGroup);
        
    }
    
    public static Map<String,ArrayList<Integer>> clusterPhaseTwo_createReadMap(ArrayList<ClusterGroup> inListGroup){
        
        /**
         * Create second map for cluster phase two purpose
         */
        
        ArrayList<Integer> dummyListNumGroup;
        Map<String,ArrayList<Integer>> readNameMap = new HashMap();
        for(int numG=0;numG<inListGroup.size();numG++){
            ClusterGroup group = inListGroup.get(numG);
            ArrayList<String> readList = group.getListReadname();
            for(int numR=0;numR<readList.size();numR++){
                String readName = readList.get(numR);
                if(readNameMap.containsKey(readName)){
                    dummyListNumGroup = readNameMap.get(readName);
                    dummyListNumGroup.add(numG);
                }else{
                    dummyListNumGroup = new ArrayList();
                    dummyListNumGroup.add(numG);
                }
                readNameMap.put(readName, dummyListNumGroup);
            }    
        }
        
        return readNameMap;
    }
    
    public static ArrayList<ClusterGroup> clusterPhaseTwo_extractCommon(ArrayList<ClusterGroup> inListGroup, Map<String,ArrayList<Integer>> inReadNameMap){
        ArrayList<ClusterGroup> listCommonGroup = new ArrayList();
        ArrayList<Integer> checkList = new ArrayList();
        for(int numG=0;numG<inListGroup.size();numG++){
            
            if(checkList.contains(numG)!=true){
                ClusterGroup group = inListGroup.get(numG);
                ArrayList<String> listRead = group.getListReadname();
                Map<Integer,ArrayList<String>> temporaryGroupMap = new HashMap();
                
                /*** Additional information ***/
                ArrayList<Byte> listNumChr = group.getListChromosome();
                ArrayList<Long> listIniPos = group.getListIniPos();
                ArrayList<Long> listLastPos = group.getListLastPos();
                ArrayList<Byte> listNumGreen = group.getListNumGreen();
                ArrayList<Byte> listNumY = group.getListNumYellow();
                ArrayList<Byte> listNumO = group.getListNumOrange();
                ArrayList<Byte> listNumR = group.getListNumRed();
                ArrayList<String> listStrand = group.getListStrand();
                ArrayList<Byte> listIniIdx = group.getListIniIndex();
                
                Map<Integer,ArrayList<Byte>> temporaryNumChrMap = new HashMap();
                Map<Integer,ArrayList<Long>> temporaryIniPosMap = new HashMap();
                Map<Integer,ArrayList<Long>> temporaryLastPosMap = new HashMap();
                Map<Integer,ArrayList<Byte>> temporaryNumGreenMap = new HashMap();
                Map<Integer,ArrayList<Byte>> temporaryNumYMap = new HashMap();
                Map<Integer,ArrayList<Byte>> temporaryNumOMap = new HashMap();
                Map<Integer,ArrayList<Byte>> temporaryNumRMap = new HashMap();
                Map<Integer,ArrayList<String>> temporaryStrandMap = new HashMap();
                Map<Integer,ArrayList<Byte>> temporaryIniIdxMap = new HashMap();
                /**********************************************************************/
                
                for(int numR=0;numR<listRead.size();numR++){
                    String dummyRead = listRead.get(numR);

                    /*** Additional information ***/
                    byte dummyNumChr = listNumChr.get(numR);
                    long dummyIniPos = listIniPos.get(numR);
                    long dummyLastPos = listLastPos.get(numR);
                    byte dummyNumGreen = listNumGreen.get(numR);
                    byte dummyNumY = listNumY.get(numR);
                    byte dummyNumO = listNumO.get(numR);        
                    byte dummyNumR = listNumR.get(numR);
                    String dummyStrand = listStrand.get(numR);
                    byte dummyIniIdx = listIniIdx.get(numR);                                
                    /**********************************************************************/
                    
                    /**
                     * Start using read name as key to get number of group from inReadNameMap
                     * and create temporary group map for each numG Group
                     */

                    ArrayList<Integer> listNumG = inReadNameMap.get(dummyRead);
                    for(int idx=0;idx<listNumG.size();idx++){
                        int dummyNumG = listNumG.get(idx);

                        if(temporaryGroupMap.containsKey(dummyNumG)){
                            ArrayList<String> dummyReadList = temporaryGroupMap.get(dummyNumG);
                            dummyReadList.add(dummyRead);
                            temporaryGroupMap.put(dummyNumG, dummyReadList);
                            
                            /*** Additional information ***/
                            ArrayList<Byte> dummyNumChrList = temporaryNumChrMap.get(dummyNumG);
                            ArrayList<Long> dummyIniPosList = temporaryIniPosMap.get(dummyNumG);
                            ArrayList<Long> dummyLastPosList = temporaryLastPosMap.get(dummyNumG);
                            ArrayList<Byte> dummyNumGreenList = temporaryNumGreenMap.get(dummyNumG);
                            ArrayList<Byte> dummyNumYList = temporaryNumYMap.get(dummyNumG);
                            ArrayList<Byte> dummyNumOList = temporaryNumOMap.get(dummyNumG);
                            ArrayList<Byte> dummyNumRList = temporaryNumRMap.get(dummyNumG);
                            ArrayList<String> dummyStrandList = temporaryStrandMap.get(dummyNumG);
                            ArrayList<Byte> dummyIniIdxList = temporaryIniIdxMap.get(dummyNumG);
                            
                            dummyNumChrList.add(dummyNumChr);
                            dummyIniPosList.add(dummyIniPos);
                            dummyLastPosList.add(dummyLastPos);
                            dummyNumGreenList.add(dummyNumGreen);
                            dummyNumYList.add(dummyNumY);
                            dummyNumOList.add(dummyNumO);
                            dummyNumRList.add(dummyNumR);
                            dummyStrandList.add(dummyStrand);
                            dummyIniIdxList.add(dummyIniIdx);
                            
                            temporaryNumChrMap.put(dummyNumG, dummyNumChrList);
                            temporaryIniPosMap.put(dummyNumG, dummyIniPosList);
                            temporaryLastPosMap.put(dummyNumG, dummyLastPosList);
                            temporaryNumGreenMap.put(dummyNumG, dummyNumGreenList);
                            temporaryNumYMap.put(dummyNumG, dummyNumYList);
                            temporaryNumOMap.put(dummyNumG, dummyNumOList);
                            temporaryNumRMap.put(dummyNumG, dummyNumRList);
                            temporaryStrandMap.put(dummyNumG, dummyStrandList);
                            temporaryIniIdxMap.put(dummyNumG, dummyIniIdxList);
                            /**********************************************************************/
                        }else{
                            ArrayList<String> dummyReadList = new ArrayList();
                            dummyReadList.add(dummyRead);
                            temporaryGroupMap.put(dummyNumG, dummyReadList);
                            
                            /*** Additional information ***/
                            ArrayList<Byte> dummyNumChrList = new ArrayList();
                            ArrayList<Long> dummyIniPosList = new ArrayList();
                            ArrayList<Long> dummyLastPosList = new ArrayList();
                            ArrayList<Byte> dummyNumGreenList = new ArrayList();
                            ArrayList<Byte> dummyNumYList = new ArrayList();
                            ArrayList<Byte> dummyNumOList = new ArrayList();
                            ArrayList<Byte> dummyNumRList = new ArrayList();
                            ArrayList<String> dummyStrandList = new ArrayList();
                            ArrayList<Byte> dummyIniIdxList = new ArrayList();
                            
                            dummyNumChrList.add(dummyNumChr);
                            dummyIniPosList.add(dummyIniPos);
                            dummyLastPosList.add(dummyLastPos);
                            dummyNumGreenList.add(dummyNumGreen);
                            dummyNumYList.add(dummyNumY);
                            dummyNumOList.add(dummyNumO);
                            dummyNumRList.add(dummyNumR);
                            dummyStrandList.add(dummyStrand);
                            dummyIniIdxList.add(dummyIniIdx);
                            
                            temporaryNumChrMap.put(dummyNumG, dummyNumChrList);
                            temporaryIniPosMap.put(dummyNumG, dummyIniPosList);
                            temporaryLastPosMap.put(dummyNumG, dummyLastPosList);
                            temporaryNumGreenMap.put(dummyNumG, dummyNumGreenList);
                            temporaryNumYMap.put(dummyNumG, dummyNumYList);
                            temporaryNumOMap.put(dummyNumG, dummyNumOList);
                            temporaryNumRMap.put(dummyNumG, dummyNumRList);
                            temporaryStrandMap.put(dummyNumG, dummyStrandList);
                            temporaryIniIdxMap.put(dummyNumG, dummyIniIdxList);
                            /**********************************************************************/
                        }
                    }
                }

                /**
                 * Find common peaks
                 */
                int oldCountRead = 0;
                ArrayList<String> selectList = new ArrayList();
                int selectNumGroup = 0;
                for(Map.Entry<Integer,ArrayList<String>> subEntry : temporaryGroupMap.entrySet()){    // minor loop for check intersect
                    Integer numGroup = subEntry.getKey();
                    ArrayList<String> listOfRead = subEntry.getValue();
                    int newCountRead = listOfRead.size();

                    if(numGroup!=numG && newCountRead>oldCountRead){
                        selectList = listOfRead;
                        oldCountRead = newCountRead;
                        
                        selectNumGroup = numGroup;
                    }
                }
            
                checkList.add(numG);
                checkList.add(selectNumGroup);
                
                if(selectList.isEmpty()!=true){
                    ClusterGroup newGroup = new ClusterGroup();
                    newGroup.addReadName(selectList);

                    /*** Additional information ***/
                    newGroup.addChromosomeNumber(temporaryNumChrMap.get(selectNumGroup));
                    newGroup.addIniPos(temporaryIniPosMap.get(selectNumGroup));
                    newGroup.addLastPos(temporaryLastPosMap.get(selectNumGroup));
                    newGroup.addNumGreen(temporaryNumGreenMap.get(selectNumGroup));
                    newGroup.addNumYellow(temporaryNumYMap.get(selectNumGroup));
                    newGroup.addNumOrange(temporaryNumOMap.get(selectNumGroup));
                    newGroup.addNumRed(temporaryNumRMap.get(selectNumGroup));
                    newGroup.addStrand(temporaryStrandMap.get(selectNumGroup));
                    newGroup.addIniIndex(temporaryIniIdxMap.get(selectNumGroup));
                    /**********************************************************************/
                    
                    listCommonGroup.add(newGroup);
                }
                
                /**** Create new group for numG group (the numG group is always one of the common peak group ****/
                ClusterGroup newGroup = new ClusterGroup();
                newGroup.addReadName(temporaryGroupMap.get(numG));

                /*** Additional information ***/
                newGroup.addChromosomeNumber(temporaryNumChrMap.get(numG));
                newGroup.addIniPos(temporaryIniPosMap.get(numG));
                newGroup.addLastPos(temporaryLastPosMap.get(numG));
                newGroup.addNumGreen(temporaryNumGreenMap.get(numG));
                newGroup.addNumYellow(temporaryNumYMap.get(numG));
                newGroup.addNumOrange(temporaryNumOMap.get(numG));
                newGroup.addNumRed(temporaryNumRMap.get(numG));
                newGroup.addStrand(temporaryStrandMap.get(numG));
                newGroup.addIniIndex(temporaryIniIdxMap.get(numG));
                /**********************************************************************/

                listCommonGroup.add(newGroup);
                /**********************************/
            }
            
        }
        return listCommonGroup;
        
    }
    
    public static ArrayList<ClusterGroup> clusterPhaseTwo_groupBySignificant(ArrayList<ClusterGroup> inListGroup){
        ArrayList<ClusterGroup> significantListGroup = new ArrayList();
        for(int numG=0;numG<inListGroup.size();numG++){
            ClusterGroup dummyGroup = inListGroup.get(numG);
            if(dummyGroup.checkSignificantFlag()==true){
                significantListGroup.add(dummyGroup);
            }
        }
        return significantListGroup;
    }
    
    
            
}
