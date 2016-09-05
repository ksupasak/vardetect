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
    private static long mask = 268435455;
    private long[] clusterVector;
    private double[] distanceVector;
    
    ArrayList<MerRead> mers;
    ArrayList listChr;                          // [long] ArrayList of align chromosome )already sort and cut with threshold [match>5 & red<5]
    ArrayList listPos;                          // [long] ArrayList of align chromosome )already sort and cut with threshold [match>5 & red<5]
    ArrayList listStrand;                       // [String] ArrayList of align chromosome )already sort and cut with threshold [match>5 & red<5]
    ArrayList listResultCode;                   // [object should cast to long (chr 5bit|strand 1bit|position 28bit)] ArrayList of align chromosome )already sort and cut with threshold [match>5 & red<5]
    ArrayList inGroup;
    ArrayList outGroup;        
    ArrayList<VectorResult> listVector;
    ArrayList<Long> frontAlgnCode;
    ArrayList<Long> backAlgnCode;
    ArrayList<ReconstructSequence> listReconSeq;
    //ArrayList<AlignmentData> algns;
    Map<Long,long[]> countResult;
    Map<Long,long[]> countResultSorted;
    Map<Long,long[]> countResultSortedCut;
    Map<Long,int[]> algnCodeIndex;             // map store the align code as key and start index , last index as value
    
    public ShortgunSequence(String seq){
        this.seq = seq;
        this.mers = new ArrayList();
        //this.algns = new ArrayList();
        this.countResult = new LinkedHashMap();
        this.countResultSorted = new LinkedHashMap();
        this.countResultSortedCut = new LinkedHashMap();
        this.algnCodeIndex = new LinkedHashMap();
        this.listChr = new ArrayList();
        this.listPos = new ArrayList();
        this.listStrand = new ArrayList();
        this.listResultCode = new ArrayList();
        this.clusterVector = new long[24];
        this.distanceVector = new double[24];
        this.inGroup = new ArrayList();
        this.outGroup = new ArrayList();
        this.frontAlgnCode = new ArrayList();
        this.backAlgnCode = new ArrayList();
        this.listReconSeq = new ArrayList();
        
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
    
    public void addDistanceVector(double[] a){
        
        /* distanceVector store the distance value of this ShortgunSequence have number of element equal to number of shortgun sequencethat has been feed in the algorithm (number of input) */
        /*  EX: input sequence have 5 shortgun sequence  */
        /*  distancevector will have 5 element also     */
        /*  this shortgun sequence is ss1                */
        /*  distanceVector = [ 0 , some distance value btw ss1 and ss2 , some distance value btw ss1 and ss3 , some distance value btw ss1 and ss4 , some distance value btw ss1 and ss5 ] */
        /* this distanceVector has been calculate from method calculateEuclidientdistancein in class AlignmentResultRead  */
        
        distanceVector = a;
    }
    
    public ArrayList<ReconstructSequence> getListReconSeq(){
        return this.listReconSeq;
    }
    
    public double[] getDistanceVector(){
        return distanceVector;
    }
    
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
    
    public ArrayList getListChrMatch(){        
        return this.listChr;
    }
    
    public ArrayList getListPosMatch(){
        return this.listPos;
    }
    
    public ArrayList getListStrand(){
        return this.listStrand;
    }
    
    public long[] getClusterVector(){ 
        /* Return vector of top two high match position of this short sequence Vector has 24 dimension equal to number of chromosome)*/
        createVectorTopTwo();
//        System.out.print(readName + ": This is clusterVector vector check:\t");
        for (int i = 0; i<clusterVector.length;i++){
//            System.out.print(clusterVector[i] + "\t");
        }
//        System.out.println();
        return clusterVector; 
    }
    
//    public void createVectorAll(String choice){
//        
//        ArrayList<VectorResult> listVector = new ArrayList();
//        for (int i=0;i<listChr.size();i++){
//            for (int j=i+1;j<listChr.size();j++){
//                VectorResult vectorResult = new VectorResult();
//                vectorResult.addChrPosCode((long)listResultCode.get(i),(long)listResultCode.get(j));
//                listVector.add(vectorResult);
//            }    
//        }
//    }
    
    public ArrayList getInGroup(){
        return this.inGroup;
    }
    
    public ArrayList getOutGroup(){
        return this.outGroup; // ArrayList outGroup is already rearrange from low to high
    }
    
    public void createInGroupOutGroup(double th){
        /* Key for clustering */
        /*                                                                                                              */
        /* Create inGroup and outGroup which is important for clustering purpose (For this shortgun Sequence)                                       */
        /* To grouping result, our assumption is: samegroup of shortgun sequence must have same outGroup                 */
        /* outGroup is ArrayList that contain the index of distance vector that have distant value above the threshold  */
        /* Ex:  we have 4 shortgunsequence (R0SS0,R0SS1,R1SS0,R1SS1)                                                    */ 
        /* R0SS0 : inGroup = 4  outGroup = 2 3                                                                          */
        /* R0SS1 : inGroup = 3  outGroup = 1 4                                                                          */
        /* R1SS0 : inGroup = 2  outGroup = 1 4                                                                          */
        /* R1SS1 : inGroup = 1  outGroup = 2 3                                                                          */
        /*                                                                                                              */
        /* You can see that we can group it to 2 group by observing the outGroup value                                   */    
        /*                                                                                                              */
        
        for(int i =0;i<this.distanceVector.length;i++){
            double dValue = this.distanceVector[i]; 
            if (dValue <= th){
                inGroup.add(i);
            }else{
                outGroup.add(i);                
            }
        }
    }
    
    public void createVectorTopTwo(){
        clusterVector = new long[24];

        listVector = new ArrayList();
        VectorResult vectorResult = new VectorResult();
        
        if(listResultCode.size()==1){
            int index = Math.toIntExact((long)listChr.get(0));
            clusterVector[index-1] = (long)listPos.get(0);          //index-1 because index in java start from 0 
        }else if(listResultCode.size()>1){                          //So for chr 24 it must be index 23 of array
            int index = Math.toIntExact((long)listChr.get(0));
            clusterVector[index-1] = (long)listPos.get(0);
            index = Math.toIntExact((long)listChr.get(1));
            clusterVector[index-1] = (long)listPos.get(1);
        }
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
        /* Contain special part that create preriquisit data for cluster propose*/
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
                    /******************************************************************/
                    /* Special part create list of chr and pos for clustering propose */
                    if(this.countResult.get(selectKey)[0]>=this.threshold && this.countResult.get(selectKey)[1]<=this.threshold){
                        this.countResultSortedCut.put((long)selectKey, this.countResult.get(selectKey));
                        long chr = ((long)selectKey)>>29;
                        long pos = ((long)selectKey&this.mask);
                        
                        String strandNot = "no";
                        if(((((long)selectKey)>>28)&1) == 1){
                            strandNot = "+";
                        }else if(((((long)selectKey)>>28)&1) == 0){
                            strandNot = "-";
                        }  
                        
                        listStrand.add(strandNot);
                        listChr.add(chr);
                        listPos.add(pos);
                        listResultCode.add(selectKey);
                        
                    }
                    /*******************************************************************/
                }

//            System.out.println("oldCount Round Loop: "+ oldCount);
            oldCount = 0;
//            System.out.println("Set zero oldCount Round Loop: "+ oldCount);
            i++;
        }
        
    }
    
    public void countAlignmentData(){
        /************** Key function to create Report ***************/
        /************************************************************/
        
//        System.out.println("");
//        System.out.println("");
        System.out.println("Do countAlignmentData on " + readName);
//        System.out.println("");
//        System.out.println("");
        
        long count;
        long green,yellow,red,orange,redInt,yellowInt,greenInt,orangeInt;
        
        for (int i=0;i<this.mers.size();i++){           // Loop Mer by Mer
//            System.out.println("\nLoop Mer by Mer :" + this.mers.get(i).getMerCode() + "Index : " + this.mers.get(i).getMerIndex());
            MerRead dummyMer = this.mers.get(i);
            /*** investigate this function ********/
            /* May cause from minus index part */
            dummyMer.createAlignmentResultStrandV2();
            /**********/
            ArrayList<Long> algnResult = dummyMer.getAlignmentResultStrand();
            
            for (int j=0;j<algnResult.size();j++){      // Loop Alignment Result
                
                long algnCode = algnResult.get(j);
                //Function get color is here
                long[] colorCode = detectColor(algnResult,j);  // A colorCode is the array of (red,yellow,orange,green,redInt,yellowInt,orangeInt,greenInt) it has three dimension
                    
//                System.out.println("This is colorCode check from function detectColor Red: "+colorCode[0]+" Yellow: "+colorCode[1]+" Orange: "+colorCode[2]+" Green: "+colorCode[3]);
//                System.out.println("This is colorCode check from function detectColor RedInt: "+colorCode[4]+" YellowInt: "+colorCode[5]+" OrangeInt: "+colorCode[6]+" GreenInt: "+colorCode[7]);

                long[] countAndColor = new long[9];
                    
                if (this.countResult.containsKey(algnCode)){

//                    System.out.println("Align At: " + algnCode);
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
//                    System.out.println("old RedInt: "+redInt);
//                    System.out.println("old yellowInt: "+yellowInt);
//                    System.out.println("old orangeInt: "+orangeInt);
//                    System.out.println("old greenInt: "+greenInt);

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
//                    System.out.println("new RedInt: "+countAndColor[5]);
//                    System.out.println("new yellowInt: "+countAndColor[6]);
//                    System.out.println("new orangeInt: "+countAndColor[7]);
//                    System.out.println("new greenInt: "+countAndColor[8]);
//
//                    System.out.println("This is countAndColor check before put to map: Align at: " + algnCode +" Count = " + count + " Red: "+countAndColor[1]+" Yellow: " + countAndColor[2] + " Orange: "+countAndColor[3]+" Green: "+ countAndColor[4]);
//                    System.out.println("This is countAndColor check before put to map: Align at: " + algnCode +" Count = " + count + " RedInt: "+countAndColor[5]+" YellowInt: " + countAndColor[6] + " OrangeInt: "+countAndColor[7]+" GreenInt: "+ countAndColor[8]);

                    this.countResult.put(algnCode, countAndColor);  
                }else{
//                    System.out.println("Do first time");
//                    System.out.println("First time Align At: " + algnCode);
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

//                    System.out.println("This is first time of countAndColor check before put to map: Align at: "+ algnCode +" Count = " + count + " Red: "+countAndColor[1]+" Yellow: " + countAndColor[2] +" Orange: "+countAndColor[3]+ " Green: "+ countAndColor[4]);
//                    System.out.println("This is first time of countAndColor check before put to map: Align at: "+ algnCode +" Count = " + count + " RedInt: "+countAndColor[5]+" YellowInt: " + countAndColor[6] +" OrangeInt: "+countAndColor[7]+ " GreenInt: "+ countAndColor[8]);

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
                long main_chrNumber = arrayCodePos.get(index)>>29; // left shift 29 bit because arrayCodePos cantain strand type bit (chr|Strand|Position)
                long compare_chrNumber = arrayCodePos.get(i)>>29;

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
    
    public void detectStrandPattern(){
        System.out.println("Detect possible pattern of " + this.readName);
        ArrayList containCheck = new ArrayList();
        
        /*  Extract beginIndex and LstIndex of each alignCode  */
        for(int i=0;i<this.mers.size();i++){                                                            // Loop Mer
            MerRead dummyMerRead = this.mers.get(i);
            ArrayList<Long> chrStrandAlgn = dummyMerRead.getAlignmentResultStrand();
            
            for(int ii=0;ii<chrStrandAlgn.size();ii++){                                                 // Loop over align code of mer
                long algnCode = chrStrandAlgn.get(ii);                                                  // get align code of mer
                
                for(int iii=0;iii<this.listResultCode.size();iii++){                                    // Loop over code cut of this shortgun Sequence
                    long algnCodeCut = (long)listResultCode.get(iii);                                   // algn code cut of shortgun sequence
                    
                    if (algnCode == algnCodeCut){                                                       // Check this mer is align or not
                        /* Match case */
                        long[] dummyCountResult = this.countResultSortedCut.get(algnCode);
                        int count = (int)dummyCountResult[0];       // get count number of this algnCode
                        
                        if (containCheck.contains(algnCode)){ 
                            
                        }else{                                                                          // First time
                            int[] Index = new int[2];
                            int beginIndex = dummyMerRead.getMerIndex();                                // calculate first and last index for this algnCode
                            int lastIndex = (beginIndex + count)-1;
                            Index[0] = beginIndex;
                            Index[1] = lastIndex;
                            System.out.println("Check beginIndex :" + Index[0]);
                            System.out.println("Check lastIndex :" + Index[1]);
                            System.out.println("algnCode " + algnCode +" chr: "+(algnCode>>29)+"Strand: " + ((algnCode>>28)&1));
                            this.algnCodeIndex.put(algnCode, Index);                                    // put algncode and Index infomation to hashmap
                            containCheck.add(algnCode);
                        }       
                    }
                }
            }
        }
        
        /*  Possible combination pattern check  */
        Set set = this.algnCodeIndex.keySet();
        Iterator mainAlgnCodeIter = set.iterator();
            
        while(mainAlgnCodeIter.hasNext()){                                                              
            Object mainAlgnCode = mainAlgnCodeIter.next();
            int mainLastIdx = 0;
            int mainBeginIdx = 0;
            int subBeginIdx = 0;
            int subLastIdx = 0;
            /* Check strand type part (in order to consider the right index)*/
            String strandMainType =null;        
            if(((((long)mainAlgnCode)>>28)&1) == 1){
                strandMainType = "+";
            }else if(((((long)mainAlgnCode)>>28)&1) == 0){
                strandMainType = "-";
            }

            if(strandMainType.equals("+")){
                mainBeginIdx = this.algnCodeIndex.get(mainAlgnCode)[0];
                mainLastIdx = this.algnCodeIndex.get(mainAlgnCode)[1];
            }else if(strandMainType.equals("-")){
                // Begin and Last must be switch
                mainLastIdx = (this.mers.size() - this.algnCodeIndex.get(mainAlgnCode)[0])-1;
                mainBeginIdx = (this.mers.size() - this.algnCodeIndex.get(mainAlgnCode)[1])-1;
            } 
                
//            int mainBeginIdx = this.algnCodeIndex.get(mainAlgnCode)[0];
//            int mainLastIdx = this.algnCodeIndex.get(mainAlgnCode)[1];
            
            Iterator subAlgnCodeIter = set.iterator();
            int checkFlag = 0;
            while(subAlgnCodeIter.hasNext()){
                Object subAlgnCode = subAlgnCodeIter.next();
//                int subBeginIdx = this.algnCodeIndex.get(subAlgnCode)[0];
//                int subLastIdx = this.algnCodeIndex.get(subAlgnCode)[1];

                /* Check strand type part (in order to consider the right index)*/
                String strandType =null;        
                if(((((long)subAlgnCode)>>28)&1) == 1){
                    strandType = "+";
                }else if(((((long)subAlgnCode)>>28)&1) == 0){
                    strandType = "-";
                }
                       
                if(strandType.equals("+")){
                    subBeginIdx = this.algnCodeIndex.get(subAlgnCode)[0];
                    subLastIdx = this.algnCodeIndex.get(subAlgnCode)[1];
                }else if(strandType.equals("-")){
                    // Begin and Last must be switch
                    subLastIdx = (this.mers.size() - this.algnCodeIndex.get(subAlgnCode)[0])-1;
                    subBeginIdx = (this.mers.size() - this.algnCodeIndex.get(subAlgnCode)[1])-1;
                } 
                /* ------------------------------ */
                System.out.println("Index Check before create new group : mainBeginIdx: "+ mainBeginIdx + "\tmainLastIdx: " + mainLastIdx + "\tsubBeginIdx: " + subBeginIdx + "\tsubLastId: " +subLastIdx);
                if (mainLastIdx < subBeginIdx){    // Case check to make sure that the front part and back part doesn't overlap                                                         // check for create new possible combination 
                    ReconstructSequence newRecon = new ReconstructSequence(this.mers,(long)mainAlgnCode,(long)subAlgnCode,mainBeginIdx,mainLastIdx,subBeginIdx,subLastIdx);                     
                    this.listReconSeq.add(newRecon);
                    checkFlag = 1;
                }
            }
            
            if(checkFlag == 0){             // checkFlag = 0 mean it's maybe single align 
                ReconstructSequence newRecon = new ReconstructSequence(this.mers,(long)mainAlgnCode,mainBeginIdx,mainLastIdx);    
                this.listReconSeq.add(newRecon);
            }
        }               
    }
                    
    public void reconstructAlignment(){
        
        ArrayList containCheck = new ArrayList();
        
        for(int i=0;i<this.mers.size();i++){                                                            // Loop Mer
            MerRead dummyMerRead = this.mers.get(i);
            ArrayList<Long> chrStrandAlgn = dummyMerRead.getAlignmentResultStrand();
            
            for(int ii=0;ii<chrStrandAlgn.size();ii++){                                                 // Loop over align code of mer
                long algnCode = chrStrandAlgn.get(ii);                                                  // get align code of mer
                
                for(int iii=0;iii<this.listResultCode.size();iii++){                                    // Loop over code cut of this shortgun Sequence
                    long algnCodeCut = (long)listResultCode.get(iii);                                   // algn code cut of shortgun sequence
                    
                    if (algnCode == algnCodeCut){                                                       // Check this mer is align or not
                        /* Match case */
                        long[] dummyCountResult = this.countResultSortedCut.get(algnCode);
                        long count = dummyCountResult[0];       // get count number of this algnCode
                        
                        dummyMerRead.addAlignSymbol("1");
                        if (containCheck.contains(algnCode)){                                           // Check for first time found (if true it is second time)
                            
                           for(int idx=0;i<listReconSeq.size();idx++){                                  // Loop in the reconstruct list to pick the right ReconstructSequence object that already be create
              
                               ReconstructSequence dummyRecon = listReconSeq.get(idx);
                               if(dummyRecon.getAlignCode().contains(algnCode)){                        // Check for the object ReconstructSequence to pick the right one (in case that they will have many object in it)
                                   
                                   ReconstructSequence selectRecon = dummyRecon;
                                   selectRecon.addMer(dummyMerRead);
                                   
                                   this.listReconSeq.add(idx, selectRecon);
                               }
                           }                                                       
                        }else{                                                                          // First time
                            /* Problen here : distingush between not contain in check because it new or not contain because it a second part */
                            ReconstructSequence selectRecon = new ReconstructSequence();
                            selectRecon.addAlignCode(algnCode);                                         // add align code to object ReconstructSequence   
                            selectRecon.addMer(dummyMerRead);                                           // add Mer
                            
                            containCheck.add(algnCode);
                            this.listReconSeq.add(selectRecon);
                        }
                        
                    }else{
                        /* non match case*/
                        dummyMerRead.addAlignSymbol("0");
                        
                        
                        
                    }
                    /* Follow note card add break point during alignment  */
                }
                
            }
        }
//        for(int ii =0; ii<this.listChr.size();ii++){
//            
//        }
//        this.listChr.

    }    
}
