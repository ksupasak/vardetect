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

}
