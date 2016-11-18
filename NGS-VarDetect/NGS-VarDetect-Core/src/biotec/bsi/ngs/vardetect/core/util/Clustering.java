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
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
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
                long mainChr = (long)matchMainChr.get(i);
                String mainStrand = (String)matchMainStrand.get(i);
                    
                for(int j=0;j<sizeMatchMainSS;j++){
                    long subChr = (long)matchSubChr.get(j);
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

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;
            String[] data = null;
            long oldIniPos = 0;
            long diff = 0;
            byte oldNumChr = 0;
            int count = 0;
            
            while ((line = reader.readLine()) != null) {
                data = line.split(",");
                
                byte numChr = Byte.parseByte(data[0]);
                long iniPos = Long.parseLong(data[1]);
//                long lastPos = Long.parseLong(data[2]);
//                byte numG = Byte.parseByte(data[3]);
//                byte numY = Byte.parseByte(data[4]);
//                byte numO = Byte.parseByte(data[5]);
//                byte numR = Byte.parseByte(data[6]);
//                String strand = data[7];
//                byte iniIdx = Byte.parseByte(data[8]);
                String readName = data[9];
                diff = iniPos - oldIniPos;
                
                if(oldNumChr == numChr){            // same Group check first criteria
                    if(diff <= maxBaseDiff){        // same Group check second criteria
                        if(group.getListReadname().contains(readName)!=true){
                            group.addReadName(readName);
                            group.addChromosomeNumber(numChr);
                            group.addIniPos(iniPos);
//                            group.addLastPos(lastPos);
//                            group.addNumGreen(numG);
//                            group.addNumYellow(numY);
//                            group.addNumOrange(numO);
//                            group.addNumRed(numR);
//                            group.addStrand(strand);
//                            group.addIniIndex(iniIdx);
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
//                        group.addLastPos(lastPos);
//                        group.addNumGreen(numG);
//                        group.addNumYellow(numY);
//                        group.addNumOrange(numO);
//                        group.addNumRed(numR);
//                        group.addStrand(strand);
//                        group.addIniIndex(iniIdx);                       
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
//                    group.addLastPos(lastPos);
//                    group.addNumGreen(numG);
//                    group.addNumYellow(numY);
//                    group.addNumOrange(numO);
//                    group.addNumRed(numR);
//                    group.addStrand(strand);
//                    group.addIniIndex(iniIdx);
                }
                
                oldNumChr = numChr;
                oldIniPos = iniPos;
                count++;
                if(count%1000000==0){
                    System.out.println(count + " line past");
                    System.out.println("Recent chromosome: " + numChr);
                }
                
            }
            writeClusterGroupToFile(filename,listGroup);
        }
        
        return listGroup;
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
            
            writer.write("Group=");
            for(int j=0; j<dummyGroup.getNumMember();j++){
//                writer.write(String.format("%d,%d,%d,%d,%d,%d,%d,%s,%d,%s", listChr.get(i),listIniPos.get(i),listLastPos.get(i),listNumG.get(i),listNumY.get(i),listNumO.get(i),listNumR.get(i),listStrand.get(i),listIniIndex.get(i),readNameList.get(i)));
                writer.write(String.format("%s",readNameList.get(i)));
                writer.write(";");
            }
            writer.write("\n");
        }
        
        writer.flush();
        writer.close();
    }
}
