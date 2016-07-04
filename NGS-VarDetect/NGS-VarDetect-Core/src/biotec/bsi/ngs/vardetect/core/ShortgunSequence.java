/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author soup
 */
public class ShortgunSequence {
 
    private String seq;
    private String readName;
    private int threshold = 0;
    private long clusterCode = 0;
    
    ArrayList<MerRead> mers;
    //ArrayList<AlignmentData> algns;
    Map<Long,long[]> countResult;
    Map<Long,long[]> countResultSorted;
    Map<Long,long[]> countResultSortedCut;
    
    public ShortgunSequence(String seq){
        this.seq = seq;
        this.mers = new ArrayList();
        //this.algns = new ArrayList();
        this.countResult = new LinkedHashMap();
        this.countResultSorted = new LinkedHashMap();
        this.countResultSortedCut = new LinkedHashMap();
    }
    
    public void addReadName(String readName){
        this.readName = readName;
    }
    
    public void addMerRead(MerRead mer){
        this.mers.add(mer);
    }
    public void addMerReadByIndex(int idx,MerRead mer){
        this.mers.set(idx, mer);
    }

//    public void addAlignmentData(AlignmentData algn){
//        this.algns.add(algn);
//    }
    
     public int getMerReadSize(){
        return this.mers.size();
    }
    
    public ArrayList<MerRead> getMerRead(){
        return this.mers;
    }
    public String getReadName(){
        return readName;
    }
    
    public String getSequence(){
        return seq;
    }
    
    public int getShortgunLength(){
        return seq.length();
    }
    
    public Map<Long,long[]> getAlignmentCount(){
        return this.countResult;
    }
    
    public Map<Long,long[]> getAlignmentCountSorted(){
        sortCountResult();
        return this.countResultSorted;
    }
    
    public Map<Long,long[]> getAlignmentCountSortedCut(int th){
        this.threshold = th;
        this.countResultSortedCut = new LinkedHashMap();
        sortCountCutResult();
        return this.countResultSortedCut;
    }
    
    public long getClusterCode(){
        this.clusterCode = 0;
        Set set = this.countResultSortedCut.keySet();
        Iterator keyIter = set.iterator();
        
        while(keyIter.hasNext()){
            //System.out.println("Read name" + readName + " Cluster Code Check : " + this.clusterCode);
            this.clusterCode = this.clusterCode + (long)keyIter.next();
            //System.out.println("Read name" + readName + " Cluster Code Check : " + this.clusterCode);
        }
        return this.clusterCode;
    }
    
    public void sortCountResult(){
        long oldCount = 0;
        long newCount = 0;
        long contaiCheck = 0;
        int i =0;
        Object selectKey = null;
        int numKey = this.countResult.size();

//        Iterator roundIter = this.countResult.keySet().iterator();
//        Iterator keyIter = this.countResult.keySet().iterator();
        
        Set set = this.countResult.keySet();
        Iterator roundIter = set.iterator();
        
//        System.out.println("Check Key Size: "+ set.size());
        
        while(roundIter.hasNext()){
              roundIter.next();
//            System.out.println("Do sorting round: "+ i);

            Iterator keyIter = set.iterator();
            while(keyIter.hasNext()){
                Object key = keyIter.next();
                if(this.countResultSorted.containsKey(key)){

                }else{
                    newCount = this.countResult.get(key)[0];
//                        System.out.println("Check newCount: "+newCount);
//                        System.out.println("Check oldCount: "+oldCount);

                    if(newCount > oldCount){
                        oldCount = newCount;
                        selectKey = key;
//                            System.out.println("New>Old (Not Exist) Check OldCount: " + oldCount);
//                            System.out.println("New>Old (Not Exist) Check Key: " + selectKey);
                    }else if(newCount == oldCount){
                        oldCount = newCount;
                        selectKey = key;
//                            System.out.println("New=Old (Not Exist) Check OldCount: " + oldCount);
//                            System.out.println("New=Old (Not Exist) Check Key: " + selectKey);
                    }
                }         
            }
//                System.out.println("Last check all value before collecting: Key="+selectKey+"Value count"+this.countResult.get(selectKey)[0]);
            if (selectKey != null){
                this.countResultSorted.put((long)selectKey, this.countResult.get(selectKey));
            }

//            System.out.println("oldCount Round Loop: "+ oldCount);
            oldCount = 0;
//            System.out.println("Set zero oldCount Round Loop: "+ oldCount);
            i++;
        }
        
    }
    
    public void sortCountCutResult(){
        long oldCount = 0;
        long newCount = 0;
        long containCheck = 0;
        int i =0;
        Object selectKey = null;
        int numKey = this.countResult.size();

//        Iterator roundIter = this.countResult.keySet().iterator();
//        Iterator keyIter = this.countResult.keySet().iterator();
        
        Set set = this.countResult.keySet();
        Iterator roundIter = set.iterator();
        
//        System.out.println("Check Key Size: "+ set.size());
        
        while(roundIter.hasNext()){
              roundIter.next();
//            System.out.println("Do sorting round: "+ i);

                Iterator keyIter = set.iterator();
                while(keyIter.hasNext()){
                    Object key = keyIter.next();
                    if(this.countResultSortedCut.containsKey(key)){

                    }else{
                        newCount = this.countResult.get(key)[0];
//                        System.out.println("Check newCount: "+newCount);
//                        System.out.println("Check oldCount: "+oldCount);

                        if(newCount > oldCount){
                            oldCount = newCount;
                            selectKey = key;
//                            System.out.println("New>Old (Not Exist) Check OldCount: " + oldCount);
//                            System.out.println("New>Old (Not Exist) Check Key: " + selectKey);
                        }else if(newCount == oldCount){
                            oldCount = newCount;
                            selectKey = key;
//                            System.out.println("New=Old (Not Exist) Check OldCount: " + oldCount);
//                            System.out.println("New=Old (Not Exist) Check Key: " + selectKey);
                        }
                    }         
                }
//                System.out.println("Last check all value before collecting: Key="+selectKey+"Value count"+this.countResult.get(selectKey)[0]);
                if(selectKey != null){
//                    if(this.countResult.get(selectKey)[0]>=this.threshold){
//                        if(this.countResult.get(selectKey)[1]<=this.threshold){
//                            
//                        }
//                        
//                    }
                    if(this.countResult.get(selectKey)[0]>=this.threshold && this.countResult.get(selectKey)[1]<=this.threshold){
                        this.countResultSortedCut.put((long)selectKey, this.countResult.get(selectKey));
                    }
                }

//            System.out.println("oldCount Round Loop: "+ oldCount);
            oldCount = 0;
//            System.out.println("Set zero oldCount Round Loop: "+ oldCount);
            i++;
        }
        
    }
    
    public void countAlignmentData(){
        System.out.println("");
        System.out.println("");
        System.out.println("Do countAlignmentData");
        System.out.println("");
        System.out.println("");
        
        long count;
        long green,yellow,red,orange,redInt,yellowInt,greenInt,orangeInt;
        
        for (int i=0;i<this.mers.size();i++){           // Loop Mer by Mer
            System.out.println("Loop Mer by Mer :" + this.mers.get(i).getMerCode() + "Index : " + this.mers.get(i).getMerIndex());
            MerRead dummyMer = this.mers.get(i);
            dummyMer.createAlignmentResult();
            ArrayList<Long> algnResult = dummyMer.getAlignmentResult();
            
            for (int j=0;j<algnResult.size();j++){      // Loop Alignment Result
                
                long algnCode = algnResult.get(j);
                //Function get color is here
                long[] colorCode = detectColor(algnResult,j);  // A colorCode is the array of (red,yellow,orange,green,redInt,yellowInt,orangeInt,greenInt) it has three dimension
                    
                System.out.println("This is colorCode check from function detectColor Red: "+colorCode[0]+" Yellow: "+colorCode[1]+" Orange: "+colorCode[2]+" Green: "+colorCode[3]);
                System.out.println("This is colorCode check from function detectColor RedInt: "+colorCode[4]+" YellowInt: "+colorCode[5]+" OrangeInt: "+colorCode[6]+" GreenInt: "+colorCode[7]);

                long[] countAndColor = new long[9];
                    
                if (this.countResult.containsKey(algnCode)){

                    System.out.println("Align At: " + algnCode);
                    countAndColor = this.countResult.get(algnCode);                       
                    count = countAndColor[0];
                    red = countAndColor[1];
                    yellow = countAndColor[2];
                    orange = countAndColor[3];
                    green = countAndColor[4];
                    redInt = countAndColor[5];
                    yellowInt = countAndColor[6];
                    orangeInt = countAndColor[7];
                    greenInt = countAndColor[8];

                    count++;                    
//                       System.out.println("old Red: "+red);
//                       System.out.println("old yellow: "+yellow);
//                       System.out.println("old orange: "+orange);
//                       System.out.println("old green: "+green);
                    System.out.println("old RedInt: "+redInt);
                    System.out.println("old yellowInt: "+yellowInt);
                    System.out.println("old orangeInt: "+orangeInt);
                    System.out.println("old greenInt: "+greenInt);

                    countAndColor[0] = count;
                    countAndColor[1] = red + colorCode[0];
                    countAndColor[2] = yellow + colorCode[1];
                    countAndColor[3] = orange + colorCode[2];
                    countAndColor[4] = green + colorCode[3];
                    countAndColor[5] = redInt + colorCode[4];
                    countAndColor[6] = yellowInt + colorCode[5];
                    countAndColor[7] = orangeInt + colorCode[6];
                    countAndColor[8] = greenInt + colorCode[7];

//                       System.out.println("new Red: "+countAndColor[1]);
//                       System.out.println("new yellow: "+countAndColor[2]);
//                       System.out.println("new orange: "+countAndColor[3]);
//                       System.out.println("new green: "+countAndColor[4]);
                    System.out.println("new RedInt: "+countAndColor[5]);
                    System.out.println("new yellowInt: "+countAndColor[6]);
                    System.out.println("new orangeInt: "+countAndColor[7]);
                    System.out.println("new greenInt: "+countAndColor[8]);

                    System.out.println("This is countAndColor check before put to map: Align at: " + algnCode +" Count = " + count + " Red: "+countAndColor[1]+" Yellow: " + countAndColor[2] + " Orange: "+countAndColor[3]+" Green: "+ countAndColor[4]);
                    System.out.println("This is countAndColor check before put to map: Align at: " + algnCode +" Count = " + count + " RedInt: "+countAndColor[5]+" YellowInt: " + countAndColor[6] + " OrangeInt: "+countAndColor[7]+" GreenInt: "+ countAndColor[8]);

                    this.countResult.put(algnCode, countAndColor);  
                }else{
                    System.out.println("Do first time");
//                       System.out.println("First time Align At: " + alignPosV2);
                    count = 1;
                    countAndColor[0] = count;
                    countAndColor[1] = colorCode[0]; // Red
                    countAndColor[2] = colorCode[1]; // Yellow
                    countAndColor[3] = colorCode[2]; // orange
                    countAndColor[4] = colorCode[3]; // green
                    countAndColor[5] = colorCode[4]; // RedInt
                    countAndColor[6] = colorCode[5]; // YellowInt
                    countAndColor[7] = colorCode[6]; // OrangeInt
                    countAndColor[8] = colorCode[7]; // GreenInt

                    System.out.println("This is first time of countAndColor check before put to map: Align at: "+ algnCode +" Count = " + count + " Red: "+countAndColor[1]+" Yellow: " + countAndColor[2] +" Orange: "+countAndColor[3]+ " Green: "+ countAndColor[4]);
                    System.out.println("This is first time of countAndColor check before put to map: Align at: "+ algnCode +" Count = " + count + " RedInt: "+countAndColor[5]+" YellowInt: " + countAndColor[6] +" OrangeInt: "+countAndColor[7]+ " GreenInt: "+ countAndColor[8]);

                    this.countResult.put(algnCode,countAndColor);
                }   
            }       
        }    
    }
        
    public long[] detectColor(ArrayList<Long> arrayCodePos, int index){
        long[] colorCode = new long[8]; //Have four color code orange = repeat in same chr; yellow = repeat with other chr; red = repeate both same and other chr; green = unique
        int arraySize = arrayCodePos.size();
        int red = 0;
        int yellow = 0;
        int orange = 0;
        int green = 0;
        int redInt = 0;
        int yellowInt = 0;
        int orangeInt = 0;
        int greenInt = 0;

        if (arraySize == 1){
            green = 1;
            greenInt++;
        }else{
            for(int i = 0; i<arraySize;i++){
                long main_chrNumber = arrayCodePos.get(index)>>28;
                long compare_chrNumber = arrayCodePos.get(i)>>28;

                if(i != index){
                    if (main_chrNumber == compare_chrNumber){
                        orange = 1;
                        orangeInt++;
                    }else if(main_chrNumber != compare_chrNumber){
                        yellow = 1;
                        yellowInt++;
                    }
                }
            }
        }
        if (orange == 1 & yellow == 1){
            red = 1;
            redInt = orangeInt+yellowInt;
            yellow = 0;
            yellowInt = 0;
            orange = 0;
            orangeInt = 0;
        }

        colorCode[0] = red;
        colorCode[1] = yellow;
        colorCode[2] = orange;
        colorCode[3] = green;
        colorCode[4] = redInt;
        colorCode[5] = yellowInt;
        colorCode[6] = orangeInt;
        colorCode[7] = greenInt;

        return colorCode;
    }
}
