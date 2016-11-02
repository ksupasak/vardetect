/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ClusterGroup;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
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
    
    public static void createColorArray(AlignmentResultRead inRes, long rLen, long merLen){
        int numMer = (int) (rLen-merLen)+1;
        int[] matchArray = new int[numMer];
        ArrayList<ShortgunSequence> listRead = inRes.getResult();
        ArrayList<int[]> listMatchArray = new ArrayList();
        Map<String,ArrayList<int[]>> colorPatternMap = new HashMap();
        
        for(int i =0;i<listRead.size();i++){ // Loop all Read
            ShortgunSequence dummySS = listRead.get(i);            
            
            ArrayList<Long> listChr = dummySS.getListChrMatch();
            ArrayList<Long> listPos = dummySS.getListPosMatch();
            ArrayList<String> listStrand = dummySS.getListStrand();
            ArrayList<Long> listNumMatch = dummySS.getListNumMatch();
            ArrayList<Long> listIniIdx = dummySS.getListIniIdx();
            long numPattern = listChr.size();
            /**
             * Start create match Array of each pattern
             *      integer array size numMer
             *  code 1 = match ; code 0 = not match
             */
            
            for(int numP=0;numP<numPattern;numP++){
                matchArray = new int[numMer];
                
                long numMatch = listNumMatch.get(numP);
                long iniIdx = listIniIdx.get(numP);
                
                for(int merIdx = (int)iniIdx;merIdx<numMatch;merIdx++){
                    matchArray[merIdx] = 1;      
                }
                listMatchArray.add(matchArray);
            }
            colorPatternMap.put(dummySS.getReadName(), listMatchArray);
        }
        
    }

}
