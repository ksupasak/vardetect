/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigInteger;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author soup
 */
public class EncodedSequence {

    //Hashtable<Long,Long> map;
    //TreeMap<Long,Long> map;
    long[] mers;            // Store reference sequence of specific chr for mapping propose
    long[] mersComp;
    long iniIndex;
    ArrayList<Long> repeatMarker;
    
    Map<Long,Long> map;
//    Map<Long,Long> repeatMarkerFront;
//    Map<Long,Long> repeatMarkerBack;
    Map<Long,Boolean> repeatIndex;
    Map<Integer,ArrayList<Integer>> linkIndexCheck;
    Map<Long,ArrayList<Integer>> alnMerMap;
    Map<Long,Long> alnCodeCheckList;
    
    long[] repeatMarkerIndex;
    int[] linkIndex;
    int mainIndex;
    int numMer;
    String name;
    String subSequence;
    boolean repeatFlag;         // use to indicate the alignment process has align with repeatMarker (repeat Found)
    
    public long mask = -268435456; // Do & operation to get mer  
    public long mask2 = 268435455; // Do & operation to get position
    public long mask36Bit = 68719476735L;
    
    
    public EncodedSequence(){
        this.mers = null;            // Store reference sequence of specific chr for mapping propose
        this.mersComp = null;
//        this.repeatMarkerFront = new LinkedHashMap();
//        this.repeatMarkerBack = new LinkedHashMap();
        this.repeatIndex = new LinkedHashMap();
        this.linkIndexCheck = new LinkedHashMap();
        this.alnMerMap = new LinkedHashMap();
        this.alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
        this.repeatFlag = false;
    }
    
    public long addStrandNotation(long in, int strand){
        
        long pos_strand = (strand<<28)+in;      // add 1 bit to indicate strand notation 0(-) or 1(+) infront of the 28 position bit
        
        return pos_strand;
    }
    
//    public void addRepeatMarkerFB(ArrayList<Map<Long,Long>> in){
//        
//        this.repeatMarkerFront = in.get(0);
//        this.repeatMarkerBack = in.get(1);
//    }
    
    public void addRepeatMarker(ArrayList<Long> in){
        
        this.repeatMarker = in;
    }
    
    public void addRepeatIndex(Map<Long,Boolean> in){
        this.repeatIndex = in;
    }
    
    public void addRepeatMarkerIndex(long[] in){
        this.repeatMarkerIndex = in;
    }
    
    public void addLinkIndex(int[] in){
        this.linkIndex = in;
    }
    
    public long[] fullAlign(long mer){  // now this function is unuse
               
        long[] listResult = align2(mer);
//        System.out.println("This is listResult check: " + listResult.length);
//        if (this.mersComp == null){
//            createComplimentStrand();
//        }
        
        long[] listResultCompliment = align2Compliment(mer);
        
        if(listResult == null && listResultCompliment != null){
            
            return listResultCompliment;
            
        }else if (listResult != null && listResultCompliment == null){  
            
            return listResult;
            
        }else if (listResult == null && listResultCompliment == null){
            
            return null;
            
        }else{
            
            int lenLR = listResult.length;
            int lenLRComp = listResultCompliment.length;
            long[] fullResult = new long[lenLR+lenLRComp];

            System.arraycopy(listResult, 0, fullResult, 0, lenLR);
            System.arraycopy(listResultCompliment, 0, fullResult, lenLR, lenLRComp);
            return fullResult;
            
        }

        
        //return listResult;
    }
    
    public long[] align2ComplimentV2(long mer){
        int strand = 0; // notation for strand -
//        createComplimentStrand(); // Caution this function will change value in mers
        int index = alignComp(mer, 0, mers.length-1);       // Call function for compliment align 
        
        int start = -1;
        int stop = -1;
        
        if(index>0){
            for(int i=index;i>=0;i--){ 
                long imer = mers[i]&mask;
                
                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }
            
            for(int i=index;i<mers.length;i++){
                long imer = mers[i]&mask;
                
                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }
            if(start<stop){
//            System.out.println(" size "+(stop-start));
                long j[] = new long[stop-start]; 
            
                for(int i =start;i<stop;i++){
                    if(i-start>=0&&i>=0)
                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
//                    System.out.println();
//                    System.out.println("Check mers the value should be 64 bit : " + mersComp[i] );
//                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
//                    System.out.println();
                }


                return j;
            }else if(start == stop){
            /**
                * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check booth value to determine this two special case.
                * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                */
               long j[] = new long[1];
               j[0] = addStrandNotation(mers[start]&mask2,strand);
               return j;

            }else{
                return null;
            } 
            
        }
        return null;
    }
    
    
    
    public long[] align2Compliment(long mer){           // this function is unuse
        int strand = 0; // notation for strand -
//        createComplimentStrand(); // Caution this function will change value in mers
        int index = alignComp(mer, 0, mersComp.length-1);       // call binary search function with initial left and right with 0 and maximum index point
        
        int start = -1;
        int stop = -1;
        
        if(index>0){
            for(int i=index;i>=0&&i>=index-200;i--){
                long imer = mersComp[i]&mask;
                
                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    
                }
            }
            
            for(int i=index;i<mersComp.length&&i<index+200;i++){
                long imer = mersComp[i]&mask;
                
                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    
                }
            }
            if(start<stop&&stop-start<500){
//            System.out.println(" size "+(stop-start));
            long j[] = new long[stop-start]; 
            
                for(int i =start;i<stop;i++){
                    if(i-start>=0&&i>=0)
                    j[i-start] = addStrandNotation(mersComp[i]&mask2,strand);
//                    System.out.println();
//                    System.out.println("Check mers the value should be 64 bit : " + mersComp[i] );
//                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
//                    System.out.println();
                }


                return j;
            }
            else
                return null;
//            System.out.println("start : "+start+" stop : "+stop+" length :"+(stop-start));
            
            
            
        }
        
        
        return null;
    }
    
    
    
    public long[] align2(long mer){
//        System.out.println("\n Do Strand + Alignment");
        int strand = 1; // Notation for strand +
        int index = align(mer, 0, mers.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
//        Old version its fix at 200 position scan up and down from middle(index)
//        int start = -1;
//        int stop = -1;
//        
//        if(index>0){
//            for(int i=index;i>=0&&i>=index-200;i--){
//                long imer = mers[i]&mask;
//                
//                if(imer!=mer){
//                    start = i+1;
//                    break;
//                }else{
//                    
//                }
//            }
//            
//            for(int i=index;i<mers.length&&i<index+200;i++){
//                long imer = mers[i]&mask;
//                
//                if(imer!=mer){
//                    stop = i;
//                    break;
//                }else{
//                    
//                }
//            }
//            if(start<stop&&stop-start<500){
////            System.out.println(" size "+(stop-start));
//            long j[] = new long[stop-start]; 
//            
//                for(int i =start;i<stop;i++){
//                    if(i-start>=0&&i>=0)
////                    j[i-start] = mers[i]&mask2;
//                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
////                    System.out.println();
////                    System.out.println("Check mers the value should be 64 bit : " + mers[i] );
////                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
////                    System.out.println();
//                }
//
//
//                return j;
//            }
//            else
//                return null;
////            System.out.println("start : "+start+" stop : "+stop+" length :"+(stop-start));
//            
//            
//            
////        }
//        
        //createComplimentStrand();
        

/**
 * New version
 */
        int start = -1;
        int stop = -1;
        
        if(index > 0){
        
            for(int i=index;i>=0;i--){ 
                long imer = mers[i]&mask;

                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }

            for(int i=index;i<mers.length;i++){
                long imer = mers[i]&mask;

                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }

            if(start<stop){
    //            System.out.println(" size "+(stop-start));
                long j[] = new long[stop-start]; 

                for(int i =start;i<stop;i++){                                       // This loop make program slow (In process to find the way to fix this)
                    if(i-start>=0&&i>=0)
                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
    //                    System.out.println();
    //                    System.out.println("Check mers the value should be 64 bit : " + mersComp[i] );
    //                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
    //                    System.out.println();
                }


                return j;
            }else if(start == stop){
                /**
                 * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                 * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check booth value to determine this two special case.
                 * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                 */
                long j[] = new long[1];
                j[0] = addStrandNotation(mers[start]&mask2,strand);
                return j;

            }else{
                return null;
            }
        }
                
        return null;
    }
    
    public long[] align3(long mer, String inSeq, int mainIdx, int numMer, Map<Integer,ArrayList<Integer>> inLinkIndexCheck){
        
        /**
         * Core function for alignment with RepeatMarker
         * This function will return long[] 64 bit compose of merCount|strand|position (less than 10 bit|1 bit|28 bit) 
         */
        int nextIndex= -1; 
        long nextMerPos = -1; 
        long nextMer = -1;
        this.subSequence = inSeq;
        this.mainIndex = mainIdx;
        this.numMer = numMer;
        Map<Integer,ArrayList<Integer>> linkIndexCheck = inLinkIndexCheck;
        if(linkIndexCheck.containsKey(mainIdx-1)){
            linkIndexCheck.remove(mainIdx-1);            // remove all linked index that coresponse to old main index
        }
        
        
        int strand = 1; // Notation for strand +
        int index = alignWithRepeatMarker(mer, 0, this.repeatMarkerIndex.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
        if(index == -1){
            return null;
        }else{
            
            /**
             * index scanning [ Scan up and down ] 
             */
            int start = -1;
            int stop = -1;
            
//            double startAon = System.currentTimeMillis();
            for(int i=index;i>=0;i--){ 
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }

            for(int i=index;i<this.repeatMarkerIndex.length;i++){
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }
//            double stopAon = System.currentTimeMillis();
//            double totalAon = stopAon - startAon;
            /***********************************************/
            
            /**
             * Small window scan and repeat continuity scan
             */
            
            if(start<stop){
    //            System.out.println(" size "+(stop-start));
//                ArrayList<Long> j = new ArrayList();                
                long j[] = new long[(stop-start)+1]; 
                int countMatch = 0;
                long merCount = 1;
                int checkRepeatNextIndex = -1;
                boolean skipFlag = false;
//                double startTime = System.currentTimeMillis();
                for(int i =start;i<=stop;i++){                                       // This loop make program slow (In process to find the way to fix this)
                    /**
                     * Do small window scan for each index
                     */
                    
                    if(linkIndexCheck.isEmpty()!=true && linkIndexCheck.containsKey(mainIdx)){
                        if(linkIndexCheck.get(mainIdx).contains(i)){
                            skipFlag = true;
                        }
                    }
                    
                    
                    if(i-start>=0 && i>=0 && skipFlag == false){
                       
                        /// Start here => repeatScan();
                        nextIndex = this.linkIndex[i];                              // i is current position that match to current mer of Big window
                        
                        for(int n=mainIdx+1;n<(inSeq.length()-numMer)+1;n++){                            // Loop for Small window scan (main index is current index from big window)
                            
                            if(nextIndex == this.mask2){
                                break;
                            }
                            nextMerPos = this.repeatMarkerIndex[nextIndex];
                            nextMer = nextMerPos&this.mask;
                            
                            String sub = inSeq.substring(n, n+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
                            long compareMer = SequenceUtil.encodeMer(sub, numMer);
                            compareMer = compareMer<<28;
                            
                            
                            
                            if(nextMer != compareMer){
                                break;
                            }else if(linkIndexCheck.isEmpty()!=true){        // check repeat of nextindex (check contain of next index in ArrayList of past next index)
                                
                                if(linkIndexCheck.containsKey(n)){
                                    if(linkIndexCheck.get(n).contains(nextIndex)){
                                       break;
                                    }else{
                                        ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                        if(dummyNextIndex==null){
                                            System.out.println("Error");
                                        }
                                        dummyNextIndex.add(nextIndex);
                                        linkIndexCheck.put(n, dummyNextIndex);
                                        
                                        nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                        merCount++;
                                    }
                                }else{
                                    
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                    
                                    nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                    merCount++;
                                }     
                            }else{
                               
                                if(linkIndexCheck.containsKey(n)){
                                    ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);  
                                }else{
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                }

                                nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                merCount++;                                         // this merCount is what we want
                            } 
                                
//                            nextMerPos = this.repeatMarkerIndex[nextIndex];
//                            nextMer = nextMerPos&this.mask;
//                            
//                            String sub = inSeq.substring(n, n+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
//                            long compareMer = SequenceUtil.encodeMer(sub, numMer);
//                            compareMer = compareMer<<28;
//                            
//                            if(nextMer != compareMer){
//                                break;
//                            }else if(this.linkIndexCheck.get(n).contains(nextIndex)){        // check repeat of nextindex (check contain of next index in ArrayList of past next index)
//                                break;
//                            }else{
//                               
//                                if(this.linkIndexCheck.containsKey(n)){
//                                    ArrayList<Integer> dummyNextIndex = this.linkIndexCheck.get(n);
//                                    dummyNextIndex.add(nextIndex);
//                                    this.linkIndexCheck.put(n, dummyNextIndex);  
//                                }else{
//                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
//                                    dummyNextIndex.add(nextIndex);
//                                    this.linkIndexCheck.put(n, dummyNextIndex);
//                                }
//
//                                nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
//                                merCount++;                                         // this merCount is what we want
//                            }
                        }
                        
//                        j[i-start] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand);
//                        j.add((merCount<<29)+addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand));
                        j[countMatch] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand);
                        countMatch++;
                    }
                    skipFlag = false;
                    merCount = 1;

                }
//                double endTime   = System.currentTimeMillis();
//                double totalTime = endTime - startTime;

                return j;
            }else if(start == stop){
                /**
                 * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                 * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check both value to determine this two special case.
                 * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                 */
                long j[] = new long[1];
//                ArrayList<Long> j = new ArrayList();
                int merCount = 1;
                j[0] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand);
//                j.add((merCount<<29)+addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand));
                return j;

            }else{
                return null;
            } 
            /*****************************************************/  
        }

    }
    
    public long[] align3Compliment(long mer, String inSeq, int mainIdx, int numMer, Map<Integer,ArrayList<Integer>> inLinkIndexCheck){
        
        /**
         * Core function for alignment with RepeatMarker
         * This function will return long[] 64 bit compose of merCount|strand|position (less than 10 bit|1 bit|28 bit) 
         */        
        int nextIndex= -1; 
        long nextMerPos = -1; 
        long nextMer = -1;
        this.subSequence = inSeq;
        this.mainIndex = mainIdx;
        this.numMer = numMer;
        Map<Integer,ArrayList<Integer>> linkIndexCheck = inLinkIndexCheck;
        if(linkIndexCheck.containsKey(mainIdx-1)){
            linkIndexCheck.remove(mainIdx-1);            // remove all linked index that coresponse to old main index
        }
        
        
        int strand = 0; // Notation for strand +
        int index = alignWithRepeatMarker(mer, 0, this.repeatMarkerIndex.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
        if(index == -1){
            return null;
        }else{
            
            /**
             * index scanning [ Scan up and down ] 
             */
            int start = -1;
            int stop = -1;

            for(int i=index;i>=0;i--){ 
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }

            for(int i=index;i<this.repeatMarkerIndex.length;i++){
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }
            /***********************************************/
            
            /**
             * Small window scan and repeat continuity scan
             */
            
            if(start<stop){
    //            System.out.println(" size "+(stop-start));
                long j[] = new long[(stop-start)+1]; 
//                ArrayList<Long> j = new ArrayList();
                int countMatch = 0;
                long merCount = 1;                                               // merCount is always start at 1 because whenever it reach this point, it's mean that at least ome mer is already match
                int checkRepeatNextIndex = -1;
                boolean skipFlag = false;
                for(int i =start;i<=stop;i++){                                       // This loop make program slow (In process to find the way to fix this)
                    /**
                     * Do small window scan for each index
                     */
                    
                    if(linkIndexCheck.isEmpty()!=true && linkIndexCheck.containsKey(mainIdx)){
                        if(linkIndexCheck.get(mainIdx).contains(i)){
                            skipFlag = true;
                        }
                    }
                    
                    if(i-start>=0 && i>=0 && skipFlag==false){
                       
                        /// Start here => repeatScan();
                        nextIndex = this.linkIndex[i];                              // i is current position that match to current mer of Big window
                        
                        for(int n=mainIdx+1;n<(inSeq.length()-numMer)+1;n++){                            // Loop for Small window scan (main index is current index from big window)

                            if(nextIndex == this.mask2){
                                break;
                            }
                            nextMerPos = this.repeatMarkerIndex[nextIndex];
                            nextMer = nextMerPos&this.mask;
                            
                            String sub = inSeq.substring(n, n+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
                            long compareMer = SequenceUtil.encodeMer(sub, numMer);
                            compareMer = compareMer<<28;
                            
                            
                            
                            if(nextMer != compareMer){
                                break;
                            }else if(linkIndexCheck.isEmpty()!=true){        // check repeat of nextindex (check contain of next index in ArrayList of past next index)
                                
                                if(linkIndexCheck.containsKey(n)){
                                    if(linkIndexCheck.get(n).contains(nextIndex)){
                                       break;
                                    }else{
                                        ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                        if(dummyNextIndex==null){
                                            System.out.println("Error");
                                        }
                                        dummyNextIndex.add(nextIndex);
                                        linkIndexCheck.put(n, dummyNextIndex);
                                        
                                        nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                        merCount++;
                                    }
                                }else{
                                    
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                    
                                    nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                    merCount++;
                                }     
                            }else{
                               
                                if(linkIndexCheck.containsKey(n)){
                                    ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);  
                                }else{
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                }

                                nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                merCount++;                                         // this merCount is what we want
                            }  
                        }
                        
//                        j[i-start] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand);
//                        j.add((merCount<<29)+addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand));
                        j[countMatch] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand);
                        countMatch++;
                    }
                    skipFlag = false;
                    merCount = 1;

                }


                return j;
            }else if(start == stop){
                /**
                 * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                 * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check both value to determine this two special case.
                 * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                 */
                long j[] = new long[1];
//                ArrayList<Long> j = new ArrayList();
                int merCount = 1;
                j[0] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand);
//                j.add((merCount<<29)+addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand));
                return j;

            }else{
                return null;
            } 
            /*****************************************************/  
        }

    }
    
    public ArrayList align4(long mer, String inSeq, int mainIdx, int numMer, Map<Integer,ArrayList<Integer>> inLinkIndexCheck,  Map<Long,Long> inAlnCodeCheckList, Map<Long,ArrayList<Integer>> inAlnMerMap){
        
        /**
         * Core function for alignment with RepeatMarker
         * This function will return long[] 64 bit compose of merCount|strand|position (less than 10 bit|1 bit|28 bit) 
         */
        ArrayList output = new ArrayList(2);
        int nextIndex= -1; 
        long nextMerPos = -1; 
        long nextMer = -1;
        this.subSequence = inSeq;
        this.mainIndex = mainIdx;
        this.numMer = numMer;
        Map<Integer,ArrayList<Integer>> linkIndexCheck = inLinkIndexCheck;
        Map<Long,ArrayList<Integer>> alnMerMap = inAlnMerMap;     // Key is align code [strand|alignposition] and value is mer code

        Map<Long,Long> alnCodeCheckList = inAlnCodeCheckList;                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
       

        if(linkIndexCheck.containsKey(mainIdx-1)){
            linkIndexCheck.remove(mainIdx-1);            // remove all linked index that coresponse to old main index
        }
        
        
        int strand = 1; // Notation for strand +
        int index = alignWithRepeatMarker(mer, 0, this.repeatMarkerIndex.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
        if(index == -1){
            this.repeatFlag = false;
            output.add(0, alnMerMap);
            output.add(1,alnCodeCheckList);
            return output;
        }else{
            
            /**
             * index scanning [ Scan up and down ] 
             */
            int start = -1;
            int stop = -1;
            
//            double startAon = System.currentTimeMillis();
            for(int i=index;i>=0;i--){ 
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }

            for(int i=index;i<this.repeatMarkerIndex.length;i++){
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }
//            double stopAon = System.currentTimeMillis();
//            double totalAon = stopAon - startAon;
            /***********************************************/
            
            /**
             * Small window scan and repeat continuity scan
             */
            
            if(start<stop){
    //            System.out.println(" size "+(stop-start));
//                ArrayList<Long> j = new ArrayList();                
                long j[] = new long[(stop-start)+1]; 
                int countMatch = 0;
                long merCount = 1;
                int checkRepeatNextIndex = -1;
                boolean skipFlag = false;
//                double startTime = System.currentTimeMillis();
                for(int i =start;i<=stop;i++){                                       // This loop make program slow (In process to find the way to fix this)
                    /**
                     * Do small window scan for each index
                     */
                    
                    if(linkIndexCheck.isEmpty()!=true && linkIndexCheck.containsKey(mainIdx)){
                        if(linkIndexCheck.get(mainIdx).contains(i)){
                            skipFlag = true;
                        }
                    }
                    
                    
                    if(i-start>=0 && i>=0 && skipFlag == false){
                       
                        /// Start here => repeatScan();
                        nextIndex = this.linkIndex[i];                              // i is current position that match to current mer of Big window
                        
                        for(int n=mainIdx+1;n<(inSeq.length()-numMer)+1;n++){                            // Loop for Small window scan (main index is current index from big window)
                            
                            if(nextIndex == this.mask2){
                                break;
                            }
                            nextMerPos = this.repeatMarkerIndex[nextIndex];
                            nextMer = nextMerPos&this.mask;
                            
                            String sub = inSeq.substring(n, n+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
                            long compareMer = SequenceUtil.encodeMer(sub, numMer);
                            compareMer = compareMer<<28;
                            
                            
                            
                            if(nextMer != compareMer){
                                break;
                            }else if(linkIndexCheck.isEmpty()!=true){        // check repeat of nextindex (check contain of next index in ArrayList of past next index)
                                
                                if(linkIndexCheck.containsKey(n)){
                                    if(linkIndexCheck.get(n).contains(nextIndex)){
                                       break;
                                    }else{
                                        ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                        if(dummyNextIndex==null){
                                            System.out.println("Error");
                                        }
                                        dummyNextIndex.add(nextIndex);
                                        linkIndexCheck.put(n, dummyNextIndex);
                                        
                                        nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                        merCount++;
                                    }
                                }else{
                                    
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                    
                                    nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                    merCount++;
                                }     
                            }else{
                               
                                if(linkIndexCheck.containsKey(n)){
                                    ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);  
                                }else{
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                }

                                nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                merCount++;                                         // this merCount is what we want
                            }        
                        }
                        
                        /**
                         * Add create alnMerMap part (We move this part from run() then put it inside align4.)
                         * (Hope it help to reduce computational time because we can get rid off one loop in run() function)
                         */
                        long mask29Bit = 536870911;
//                        int merCount = (int)(posR[j]>>29);       
                        long alnCode = addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand) - index;     // posR is ~39 bit [merCount|strand|position] ; algncode is 29 bit [strand|alignPosition]. alignposition is position - index


                        if(alnCodeCheckList.containsKey(alnCode)){

                            iniIndex = alnCodeCheckList.get(alnCode);

                            long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                            ArrayList<Integer> merList = alnMerMap.get(indexAlnCode);

                            /**
                             * Case check to solve the problem. In case, when position-index is the same value but actually it different peak.
                             * To check continuity of this alnCode. We reserve index 0 of merList to store the recent index.
                             * Check continuity of index from different between recent index and current index.
                             */
                            if(merList==null){
                                System.out.println();
                            }
                            if(index-merList.get(0)==1){                                // Case check to solve the problem. In case, when position-index is the same value but actually it different peak
                                /**
                                 * it's continue. So, iniIndex not change 
                                 */

                                iniIndex = alnCodeCheckList.get(alnCode);

                                indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                merList.remove(0);
                                merList.add(0,index);
                                for(int num=0;num<merCount;num++){
                                    merList.add(1);
                                }

                                alnMerMap.put(indexAlnCode, merList);
                            }else{
                                /**
                                 * it's not continue. So, iniIndex has change to present index                                                                                  
                                 */

                                iniIndex = index;
                                alnCodeCheckList.put(alnCode, iniIndex);

                                indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                merList = new ArrayList();
                                merList.add(0,index);
                                for(int num=0;num<merCount;num++){
                                    merList.add(1);
                                }
                                alnMerMap.put(indexAlnCode,merList);
                            }
                            /**************************************************************************************************/

                        }else{
                            iniIndex = index;
                            alnCodeCheckList.put(alnCode, iniIndex);

                            long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                            ArrayList<Integer> merList = new ArrayList();
                            merList.add(0,index);
                            for(int num=0;num<merCount;num++){
                                merList.add(1);
                            }
                            alnMerMap.put(indexAlnCode,merList);
                        }
                        
//                        j[countMatch] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand);
//                        countMatch++;
                    }
                    skipFlag = false;
                    merCount = 1;

                }
//                double endTime   = System.currentTimeMillis();
//                double totalTime = endTime - startTime;
                this.repeatFlag = true;
                output.add(0, alnMerMap);
                output.add(1,alnCodeCheckList);
                return output; 
            }else if(start == stop){
                /**
                 * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                 * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check both value to determine this two special case.
                 * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                 */
                long j[] = new long[1];
//                ArrayList<Long> j = new ArrayList();
                int merCount = 1;
                
                long mask29Bit = 536870911;
//                        int merCount = (int)(posR[j]>>29);       
                long alnCode = addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand) - index;     // posR is ~39 bit [merCount|strand|position] ; algncode is 29 bit [strand|alignPosition]. alignposition is position - index
                
                if(alnCodeCheckList.containsKey(alnCode)){

                    iniIndex = alnCodeCheckList.get(alnCode);

                    long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                    ArrayList<Integer> merList = alnMerMap.get(indexAlnCode);

                    /**
                     * Case check to solve the problem. In case, when position-index is the same value but actually it different peak.
                     * To check continuity of this alnCode. We reserve index 0 of merList to store the recent index.
                     * Check continuity of index from different between recent index and current index.
                     */

                    if(index-merList.get(0)==1){                                // Case check to solve the problem. In case, when position-index is the same value but actually it different peak
                        /**
                         * it's continue. So, iniIndex not change 
                         */

                        iniIndex = alnCodeCheckList.get(alnCode);

                        indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                        merList.remove(0);
                        merList.add(0,index);
                        for(int num=0;num<merCount;num++){
                            merList.add(1);
                        }

                        alnMerMap.put(indexAlnCode, merList);
                    }else{
                        /**
                         * it's not continue. So, iniIndex has change to present index                                                                                  
                         */

                        iniIndex = index;
                        alnCodeCheckList.put(alnCode, iniIndex);

                        indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                        merList = new ArrayList();
                        merList.add(0,index);
                        for(int num=0;num<merCount;num++){
                            merList.add(1);
                        }
                        alnMerMap.put(indexAlnCode,merList);
                    }
                    /**************************************************************************************************/

                }else{
                    iniIndex = index;
                    alnCodeCheckList.put(alnCode, iniIndex);

                    long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                    ArrayList<Integer> merList = new ArrayList();
                    merList.add(0,index);
                    for(int num=0;num<merCount;num++){
                        merList.add(1);
                    }
                    alnMerMap.put(indexAlnCode,merList);
                }
                
                
                
//                j[0] = (merCount<<29)+addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand);
//                j.add((merCount<<29)+addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand));
                this.repeatFlag = true;
                output.add(0, alnMerMap);
                output.add(1,alnCodeCheckList);
                return output;

            }else{
                this.repeatFlag = false;
                output.add(0, alnMerMap);
                output.add(1,alnCodeCheckList);
                return output;
            } 
            /*****************************************************/  
        }

    }
    
    public ArrayList align4Compliment(long mer, String inSeq, int mainIdx, int numMer, Map<Integer,ArrayList<Integer>> inLinkIndexCheck, Map<Long,Long> inAlnCodeCheckList, Map<Long,ArrayList<Integer>> inAlnMerMap){
        
        /**
         * Core function for alignment with RepeatMarker
         * This function will return long[] 64 bit compose of merCount|strand|position (less than 10 bit|1 bit|28 bit) 
         */
        ArrayList output = new ArrayList(2);        
        int nextIndex= -1; 
        long nextMerPos = -1; 
        long nextMer = -1;
        this.subSequence = inSeq;
        this.mainIndex = mainIdx;
        this.numMer = numMer;
        Map<Integer,ArrayList<Integer>> linkIndexCheck = inLinkIndexCheck;
        Map<Long,ArrayList<Integer>> alnMerMap = inAlnMerMap;     // Key is align code [strand|alignposition] and value is mer code

        Map<Long,Long> alnCodeCheckList = inAlnCodeCheckList;                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
        
        if(linkIndexCheck.containsKey(mainIdx-1)){
            linkIndexCheck.remove(mainIdx-1);            // remove all linked index that coresponse to old main index
        }
        
        
        int strand = 0; // Notation for strand +
        int index = alignWithRepeatMarker(mer, 0, this.repeatMarkerIndex.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
        if(index == -1){
            this.repeatFlag = false;
            output.add(0, alnMerMap);
            output.add(1,alnCodeCheckList);
            return output;
        }else{
            
            /**
             * index scanning [ Scan up and down ] 
             */
            int start = -1;
            int stop = -1;

            for(int i=index;i>=0;i--){ 
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }

            for(int i=index;i<this.repeatMarkerIndex.length;i++){
                long imer = this.repeatMarkerIndex[i]&mask;

                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }
            /***********************************************/
            
            /**
             * Small window scan and repeat continuity scan
             */
            
            if(start<stop){
    //            System.out.println(" size "+(stop-start));
                long j[] = new long[(stop-start)+1]; 
//                ArrayList<Long> j = new ArrayList();
                int countMatch = 0;
                long merCount = 1;                                               // merCount is always start at 1 because whenever it reach this point, it's mean that at least ome mer is already match
                int checkRepeatNextIndex = -1;
                boolean skipFlag = false;
                for(int i =start;i<=stop;i++){                                       // This loop make program slow (In process to find the way to fix this)
                    /**
                     * Do small window scan for each index
                     */
                    
                    if(linkIndexCheck.isEmpty()!=true && linkIndexCheck.containsKey(mainIdx)){
                        if(linkIndexCheck.get(mainIdx).contains(i)){
                            skipFlag = true;
                        }
                    }
                    
                    if(i-start>=0 && i>=0 && skipFlag==false){
                       
                        /// Start here => repeatScan();
                        nextIndex = this.linkIndex[i];                              // i is current position that match to current mer of Big window
                        
                        for(int n=mainIdx+1;n<(inSeq.length()-numMer)+1;n++){                            // Loop for Small window scan (main index is current index from big window)

                            if(nextIndex == this.mask2){
                                break;
                            }
                            nextMerPos = this.repeatMarkerIndex[nextIndex];
                            nextMer = nextMerPos&this.mask;
                            
                            String sub = inSeq.substring(n, n+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
                            long compareMer = SequenceUtil.encodeMer(sub, numMer);
                            compareMer = compareMer<<28;
                            
                            
                            
                            if(nextMer != compareMer){
                                break;
                            }else if(linkIndexCheck.isEmpty()!=true){        // check repeat of nextindex (check contain of next index in ArrayList of past next index)
                                
                                if(linkIndexCheck.containsKey(n)){
                                    if(linkIndexCheck.get(n).contains(nextIndex)){
                                       break;
                                    }else{
                                        ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                        if(dummyNextIndex==null){
                                            System.out.println("Error");
                                        }
                                        dummyNextIndex.add(nextIndex);
                                        linkIndexCheck.put(n, dummyNextIndex);
                                        
                                        nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                        merCount++;
                                    }
                                }else{
                                    
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                    
                                    nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                    merCount++;
                                }     
                            }else{
                               
                                if(linkIndexCheck.containsKey(n)){
                                    ArrayList<Integer> dummyNextIndex = linkIndexCheck.get(n);
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);  
                                }else{
                                    ArrayList<Integer> dummyNextIndex = new ArrayList();
                                    dummyNextIndex.add(nextIndex);
                                    linkIndexCheck.put(n, dummyNextIndex);
                                }

                                nextIndex = this.linkIndex[nextIndex];              // update next index with old nextIndex
                                merCount++;                                         // this merCount is what we want
                            }  
                        }
                        
                        /**
                         *  Add create alnMerMap
                         */
                        long alnCode = addStrandNotation(this.repeatMarkerIndex[i]&mask2,strand) - index;     // posR is ~39 bit [merCount|strand|position] ; algncode is 29 bit [strand|alignPosition]. alignposition is position - index


                        if(alnCodeCheckList.containsKey(alnCode)){

                            iniIndex = alnCodeCheckList.get(alnCode);

                            long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                            ArrayList<Integer> merList = alnMerMap.get(indexAlnCode);

                            /**
                             * Case check to solve the problem. In case, when position-index is the same value but actually it different peak.
                             * To check continuity of this alnCode. We reserve index 0 of merList to store the recent index.
                             * Check continuity of index from different between recent index and current index.
                             */

                            if(index-merList.get(0)==1){                                // Case check to solve the problem. In case, when position-index is the same value but actually it different peak
                                /**
                                 * it's continue. So, iniIndex not change 
                                 */

                                iniIndex = alnCodeCheckList.get(alnCode);

                                indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                merList.remove(0);
                                merList.add(0,index);
                                for(int num=0;num<merCount;num++){
                                    merList.add(1);
                                }

                                alnMerMap.put(indexAlnCode, merList);
                            }else{
                                /**
                                 * it's not continue. So, iniIndex has change to present index                                                                                  
                                 */

                                iniIndex = index;
                                alnCodeCheckList.put(alnCode, iniIndex);

                                indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                merList = new ArrayList();
                                merList.add(0,index);
                                for(int num=0;num<merCount;num++){
                                    merList.add(1);
                                }
                                alnMerMap.put(indexAlnCode,merList);
                            }
                            /**************************************************************************************************/

                        }else{
                            iniIndex = index;
                            alnCodeCheckList.put(alnCode, iniIndex);

                            long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                            ArrayList<Integer> merList = new ArrayList();
                            merList.add(0,index);
                            for(int num=0;num<merCount;num++){
                                merList.add(1);
                            }
                            alnMerMap.put(indexAlnCode,merList);

                        }
                        
                    }
                    skipFlag = false;
                    merCount = 1;

                }

                this.repeatFlag = true;
                output.add(0, alnMerMap);
                output.add(1,alnCodeCheckList);
                return output;
            }else if(start == stop){
                /**
                 * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                 * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check both value to determine this two special case.
                 * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                 */
                long j[] = new long[1];
//                ArrayList<Long> j = new ArrayList();
                int merCount = 1;
                           
                long alnCode = addStrandNotation(this.repeatMarkerIndex[start]&mask2,strand) - index;     // posR is ~39 bit [merCount|strand|position] ; algncode is 29 bit [strand|alignPosition]. alignposition is position - index


                if(alnCodeCheckList.containsKey(alnCode)){

                    iniIndex = alnCodeCheckList.get(alnCode);

                    long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                    ArrayList<Integer> merList = alnMerMap.get(indexAlnCode);

                    /**
                     * Case check to solve the problem. In case, when position-index is the same value but actually it different peak.
                     * To check continuity of this alnCode. We reserve index 0 of merList to store the recent index.
                     * Check continuity of index from different between recent index and current index.
                     */

                    if(index-merList.get(0)==1){                                // Case check to solve the problem. In case, when position-index is the same value but actually it different peak
                        /**
                         * it's continue. So, iniIndex not change 
                         */

                        iniIndex = alnCodeCheckList.get(alnCode);

                        indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                        merList.remove(0);
                        merList.add(0,index);
                        for(int num=0;num<merCount;num++){
                            merList.add(1);
                        }

                        alnMerMap.put(indexAlnCode, merList);
                    }else{
                        /**
                         * it's not continue. So, iniIndex has change to present index                                                                                  
                         */

                        iniIndex = index;
                        alnCodeCheckList.put(alnCode, iniIndex);

                        indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                        merList = new ArrayList();
                        merList.add(0,index);
                        for(int num=0;num<merCount;num++){
                            merList.add(1);
                        }
                        alnMerMap.put(indexAlnCode,merList);
                    }
                    /**************************************************************************************************/

                }else{
                    iniIndex = index;
                    alnCodeCheckList.put(alnCode, iniIndex);

                    long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                    ArrayList<Integer> merList = new ArrayList();
                    merList.add(0,index);
                    for(int num=0;num<merCount;num++){
                        merList.add(1);
                    }
                    alnMerMap.put(indexAlnCode,merList);

                }
                
                this.repeatFlag = true;
                output.add(0, alnMerMap);
                output.add(1,alnCodeCheckList);
                return output;

            }else{
                this.repeatFlag = false;
                output.add(0, alnMerMap);
                output.add(1,alnCodeCheckList);
                return output;
            } 
            /*****************************************************/  
        }

    }
    
    public long[] align5(long mer){
//        System.out.println("\n Do Strand + Alignment");
        int strand = 1; // Notation for strand +
        int index = align(mer, 0, mers.length-1); // call binary search function with initial left and right with 0 and maximum index point

        /**
         * New version
         */
        int start = -1;
        int stop = -1;
        
        if(index > 0){
        
            for(int i=index;i>=0&&i>=index-1;i--){ 
                long imer = mers[i]&mask;

                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }

            for(int i=index;i<mers.length&&i<index+1;i++){
                long imer = mers[i]&mask;

                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }

            if(start<stop){
    //            System.out.println(" size "+(stop-start));
                long j[] = new long[stop-start]; 

                for(int i =start;i<stop;i++){                                       // This loop make program slow (In process to find the way to fix this)
                    if(i-start>=0&&i>=0)
                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
    //                    System.out.println();
    //                    System.out.println("Check mers the value should be 64 bit : " + mersComp[i] );
    //                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
    //                    System.out.println();
                }


                return j;
            }else if(start == stop){
                /**
                 * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                 * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check booth value to determine this two special case.
                 * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                 */
                long j[] = new long[1];
                j[0] = addStrandNotation(mers[start]&mask2,strand);
                return j;

            }else{
                return null;
            }
        }
                
        return null;
    }
    
    public long[] align5Compliment(long mer){
        int strand = 0; // notation for strand -
//        createComplimentStrand(); // Caution this function will change value in mers
        int index = alignComp(mer, 0, mers.length-1);       // Call function for compliment align 
        
        int start = -1;
        int stop = -1;
        
        if(index>0){
            for(int i=index;i>=0&&i>=index-1;i--){ 
                long imer = mers[i]&mask;
                
                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    start = i;
                }
            }
            
            for(int i=index;i<mers.length&&i<index+200;i++){
                long imer = mers[i]&mask;
                
                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    stop = i;
                }
            }
            if(start<stop){
//            System.out.println(" size "+(stop-start));
                long j[] = new long[stop-start]; 
            
                for(int i =start;i<stop;i++){
                    if(i-start>=0&&i>=0)
                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
//                    System.out.println();
//                    System.out.println("Check mers the value should be 64 bit : " + mersComp[i] );
//                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
//                    System.out.println();
                }


                return j;
            }else if(start == stop){
            /**
                * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
                * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check booth value to determine this two special case.
                * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
                */
               long j[] = new long[1];
               j[0] = addStrandNotation(mers[start]&mask2,strand);
               return j;

            }else{
                return null;
            } 
            
        }
        return null;
    }
    
    public long[] align6(long mer){
        
        /**
         * It take a roll to search mer in repeatMarker. Just check the mer is repeat or not
         * This function will return null if it not repeat and return -1 if it repeat 
         */
//        int nextIndex= -1; 
//        long nextMerPos = -1; 
//        long nextMer = -1;
//        this.subSequence = inSeq;
//        this.mainIndex = mainIdx;
//        this.numMer = numMer;
//        Map<Integer,ArrayList<Integer>> linkIndexCheck = inLinkIndexCheck;
//        if(linkIndexCheck.containsKey(mainIdx-1)){
//            linkIndexCheck.remove(mainIdx-1);            // remove all linked index that coresponse to old main index
//        }
        
        
//        int strand = 1; // Notation for strand +
        int index = alignWithRepeatMarker(mer, 0, this.repeatMarkerIndex.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
        if(index == -1){
            return null;
        }else{
            long[] j = new long[1];
            j[0] = -1;
            return j;
        }

    }
    
    public long[] align6Compliment(long mer){
        
        /**
         * It take a roll to search mer in repeatMarker. Just check the mer is repeat or not
         * This function will return null if it not repeat and return -1 if it repeat 
         */        
//        int nextIndex= -1; 
//        long nextMerPos = -1; 
//        long nextMer = -1;
//        this.subSequence = inSeq;
//        this.mainIndex = mainIdx;
//        this.numMer = numMer;
//        Map<Integer,ArrayList<Integer>> linkIndexCheck = inLinkIndexCheck;
//        if(linkIndexCheck.containsKey(mainIdx-1)){
//            linkIndexCheck.remove(mainIdx-1);            // remove all linked index that coresponse to old main index
//        }
        
        
//        int strand = 0; // Notation for strand +
        int index = alignWithRepeatMarker(mer, 0, this.repeatMarkerIndex.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
        if(index == -1){
            return null;
        }else{
            
            long[] j = new long[1];
            j[0] = -1;
            return j;
            
        }

    }
    
    public long[] alignFullMerPos(long mer){
        /**
         * this function will return index on mers array (Reference chr.bin)
         * Has been use for create repeat marker
         */
         
        
           
        
//        System.out.println("\n Do Strand + Alignment");
        int strand = 1; // Notation for strand +
        int index = alignMerPos(mer, 0, mers.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
//        Old version its fix at 200 position scan up and down from middle(index)
//        int start = -1;
//        int stop = -1;
//        
//        if(index>0){
//            for(int i=index;i>=0&&i>=index-200;i--){
//                long imer = mers[i]&mask;
//                
//                if(imer!=mer){
//                    start = i+1;
//                    break;
//                }else{
//                    
//                }
//            }
//            
//            for(int i=index;i<mers.length&&i<index+200;i++){
//                long imer = mers[i]&mask;
//                
//                if(imer!=mer){
//                    stop = i;
//                    break;
//                }else{
//                    
//                }
//            }
//            if(start<stop&&stop-start<500){
////            System.out.println(" size "+(stop-start));
//            long j[] = new long[stop-start]; 
//            
//                for(int i =start;i<stop;i++){
//                    if(i-start>=0&&i>=0)
////                    j[i-start] = mers[i]&mask2;
//                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
////                    System.out.println();
////                    System.out.println("Check mers the value should be 64 bit : " + mers[i] );
////                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
////                    System.out.println();
//                }
//
//
//                return j;
//            }
//            else
//                return null;
////            System.out.println("start : "+start+" stop : "+stop+" length :"+(stop-start));
//            
//            
//            
////        }
//        
        //createComplimentStrand();
        

/**
 * New version
 */
        int start = -1;
        int stop = -1;
        
                               
        for(int i=index;i>=0;i--){                   // Scan up
            long imer = mers[i];

            if(imer!=mer){
                start = i+1;
                break;
            }else{
                start = i;
            }
        }

        for(int i=index;i<mers.length;i++){     // Scan down
            long imer = mers[i];

            if(imer!=mer){
                stop = i;
                break;
            }else{
                stop = i;
            }
        }
        
        if(start<stop){
//            System.out.println(" size "+(stop-start));
            long j[] = new long[stop-start]; 

            for(int i =start;i<stop;i++){
                if(i-start>=0&&i>=0)
                j[i-start] = i;

            }


            return j;
        }else if(start == stop){
            /**
             * In case of index has value equal to 0 and 28bit value. We cannot scan up for the case that index is 0 and we cannot scan down for the case that index is 28bit(maximum index)
             * With this two case the scan protocol above will return the same value of start and stop index If it not repeat. So, we can check booth value to determine this two special case.
             * If it repeat it will fall into above check case (Because, with the repeat we can possibly scan down or scan up).
             */
            long j[] = new long[1];
            j[0] = addStrandNotation(mers[start]&mask2,strand);
            return j;
        }else{
            return null; 
        }
               
//            System.out.println("start : "+start+" stop : "+stop+" length :"+(stop-start));

    }
 
    public long[] alignLocal(long mer){
//        System.out.println("\n Do Strand + Alignment");
        int strand = 1; // Notation for strand +
        int index = align(mer, 0, mers.length-1); // call binary search function with initial left and right with 0 and maximum index point
        
        int start = -1;
        int stop = -1;
        
        if(index>0){
            for(int i=index;i>=0&&i>=index-200;i--){
                long imer = mers[i]&mask;
                
                if(imer!=mer){
                    start = i+1;
                    break;
                }else{
                    
                }
            }
            
            for(int i=index;i<mers.length&&i<index+200;i++){
                long imer = mers[i]&mask;
                
                if(imer!=mer){
                    stop = i;
                    break;
                }else{
                    
                }
            }
            if(start<stop&&stop-start<500){
//            System.out.println(" size "+(stop-start));
            long j[] = new long[stop-start]; 
            
                for(int i =start;i<stop;i++){
                    if(i-start>=0&&i>=0)
//                    j[i-start] = mers[i]&mask2;
                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
//                    System.out.println();
//                    System.out.println("Check mers the value should be 64 bit : " + mers[i] );
//                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
//                    System.out.println();
                }


                return j;
            }
            else
                return null;
//            System.out.println("start : "+start+" stop : "+stop+" length :"+(stop-start));
            
            
            
        }
        
        //createComplimentStrand();
        return null;
    }
    
    public long align(long mer){        // Curently unuse
        
        int index = align(mer, 0, mers.length-1);
        if(index>0)return mers[index]&mask2;
        return -1;
    }
    
    public int align(long mer, int left, int right){
        // Core of Binary Search Function (normal strand(+))
        int mid = (left+right)/2;       // Find middle point between left and right
        long i = mers[mid]&mask;        // get reference mers code at that middle point for matching purpose  !! Important reference mer must be sorted 
        
        if(left>right)return -1;        // in case that left value higher than right that mean this mer not match
        else
            if(i<mer){                  // if selected reference mers code less than input mer
                return align(mer, mid+1,right);     // adjust to new index by subtitude left position to mid+1 
            }else
                if(i>mer){              // if selected reference mers code higher than input mer
                    return align(mer, left,mid-1);  // adjust to new index by subtitude right position to mid-1
                }else
                    if(i==mer){         // if equal mean this position is match
//                        long j = mers[mid];
                        return mid;     // return match position
                    }
        
        return -1;
    }
    
    public int alignWithRepeatMarker(long mer, int left, int right){
        // Core of Binary Search Function (normal strand(+))
        long dummyMer = mer&mask;
        int mid = (left+right)/2;       // Find middle point between left and right
        long i = this.repeatMarkerIndex[mid]&mask;        // get reference mers code at that middle point for matching purpose  !! Important reference mer must be sorted 

        if(left>right)return -1;        // in case that left value higher than right that mean this mer not match
        else
            if(i<mer){                  // if selected reference mers code less than input mer
                return alignWithRepeatMarker(mer, mid+1,right);     // adjust to new index by subtitude left position to mid+1 
            }else
                if(i>mer){              // if selected reference mers code higher than input mer
                    return alignWithRepeatMarker(mer, left,mid-1);  // adjust to new index by subtitude right position to mid-1
                }else
                    if(i==mer){         // if equal mean this position is match
    //                        long j = mers[mid];
                        return mid;     // return match position
                    }       
        return -1;
    }
    
    public int alignMerPos(long mer, int left, int right){
        // Core of Binary Search Function (normal strand(+))
        int mid = (left+right)/2;       // Find middle point between left and right
        long i = mers[mid];        // get reference mers code at that middle point for matching purpose (long i = mer|Pos)  !! Important reference mer must be sorted 
        
        if(left>right)return -1;        // in case that left value higher than right that mean this mer not match
        else
            if(i<mer){                  // if selected reference mers code less than input mer
                return alignMerPos(mer, mid+1,right);     // adjust to new index by subtitude left position to mid+1 
            }else
                if(i>mer){              // if selected reference mers code higher than input mer
                    return alignMerPos(mer, left,mid-1);  // adjust to new index by subtitude right position to mid-1
                }else
                    if(i==mer){         // if equal mean this position is match
//                        long j = mers[mid];
                        return mid;     // return match position
                    }
        
        return -1;
    }
    
    public long alignComp(long mer){        // Currently not use
        
        int index = alignComp(mer, 0, mers.length-1);
        if(index>0)return mers[index]&mask2;
        return -1;
    }
    
    public int alignComp(long mer, int left, int right){
        // Core of Binary Search Function (compliment strand(-))
        int mid = (left+right)/2;
        long i = mers[mid]&mask;
        
        if(left>right)return -1;
        else
            if(i<mer){
                return alignComp(mer, mid+1,right);
            }else
                if(i>mer){
                    return alignComp(mer, left,mid-1);
                }else
                    if(i==mer){
//                        long j = mers[mid];
                        return mid;
                    }
        
        return -1;
    }
    
    public void setMersComp(long mers[]){       // Curently not use mersComp (Use on old Implementation)
        this.mersComp = mers;     
    }
    
    public void setMers(long mers[]){
        this.mers = mers;     
    }
    
    public long [] getMersComp(){
        return this.mersComp;
    }
    
    public long [] getMers(){
        return this.mers;
    }
    
    public boolean getRepeatFlag(){
        return this.repeatFlag;
    }
    
    public void setMap(Map<Long, Long> map) {       //Currently not use (Use on old implementation hashmap version)
        this.map = map;
        //this.name = chrName;
    }
    
    public void setReadMap(Map<Long,Long> map){     //Currently not use (Use on old implementation hashmap version)
        this.map = map;
    }
    
    public Map getEncodeMap(){                      //Currently not use (Use on old implementation hashmap version)
        return this.map;
    }
    
    public String getEncodeChrName(){               //Currently not use 
        return this.name;
    }
    
    public void readFromPath(String file_path, String fa) throws FileNotFoundException, IOException {       // Currently not use (use in old implementation)
        
        map = new HashMap<Long,Long>();
        
        
        if(fa.compareTo("map")==0){
        
            Charset charset = Charset.forName("US-ASCII");

            Path path = Paths.get(file_path+"."+fa);


            StringBuffer seq = new StringBuffer();

            try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
                String line = null;
                int count = 0 ;



                while ((line = reader.readLine()) != null) {
                    String[] st = line.split("\t");

                    long mer = Long.valueOf(st[0]);
                    long pos = Long.valueOf(st[1]);

                    map.put(mer, pos);

                    if(count%1000000==0)System.out.println("Read Mer "+count);
                    count ++;
                }
            System.out.println("Total mer : "+map.size());
        
            }catch(Exception e){
           
            }
        
        }else if(fa.compareTo("bmap")==0){
            int count = 0 ;
            
            DataInputStream is = new DataInputStream(new FileInputStream(file_path+"."+fa));
            int size = is.readInt();
            System.out.println("Totalxx bmer : "+size);

            for(int i=0;i<size;i++){
               
                long mer = is.readLong();
                long pos = is.readLong();
                map.put(mer, pos);
            
                if(count%1000000==0)System.out.println("Read binary Mer "+count);
                    count ++;
            }
            System.out.println("Total bmer : "+size);

            is.close();
        }else if(fa.compareTo("bin")==0){
            int count = 0 ;
            DataInputStream is = new DataInputStream(new FileInputStream(file_path+"."+fa));
            
            int size = is.readInt();
            System.out.println("Totalxx bmer : "+size);
            //System.out.println(is.readLong()>>28);
            //System.out.println(is.readLong()&268435455);

            for(int i=0;i<size;i++){
                
                long mer = is.readLong();
                //mers[i] = mer;
                
                if((mer&268435455)<100){
                     System.out.println("Yeahhhhhhhh" + (mer&268435455));
                }
                if(count%1000000==0)System.out.println("Read binary Mer "+count);
                    count ++;
            }
            System.out.println("Total bmer : "+size);

            is.close();  
        }
      
        
    }
    
    
    public void writeToPath(String path, String fa) throws FileNotFoundException, IOException {         //Currently not use (Use on old implementation hashmap version)

       
       
        //Enumeration<Long> e = map.keys();
        if(fa.compareTo("map")==0){
            PrintStream ps = new PrintStream(path+"."+fa);
            for (Map.Entry<Long,Long> entry : map.entrySet()){
                Long mer = entry.getKey();
                Long pos = map.get(mer);
                ps.println(mer+"\t"+pos);

            }
        }
        else if(fa.compareTo("bmap")==0){

            DataOutputStream os = new DataOutputStream(new FileOutputStream(path+"."+fa));
            System.out.println("Total bmer : "+map.keySet().size());

            os.writeInt(map.keySet().size());
            for (Map.Entry<Long,Long> entry : map.entrySet()){
                Long mer = entry.getKey();
                Long pos = map.get(mer);
                os.writeLong(mer);
                os.writeLong(pos);
            }
            os.close();    
        } 

    }
       
       
       
       /*while(e.hasMoreElements()){
           Long mer = e.nextElement();
           Long pos = map.get(mer);
           ps.println(mer+"\t"+pos);
         
       }*/

    public void lazyLoad() {
        
         this.mers = null;
         this.mersComp = null;

    }
    
    public void createComplimentStrand(){           // Currently not use (use on do complement reference version)
        
        System.out.println("\n Create compliment strand ");
        this.mersComp = Arrays.copyOf(mers, mers.length);
        //this.mersComp = this.mers;
        for(int i=0;i<this.mersComp.length;i++){
            long dummyMerPos = this.mersComp[i];
//            System.out.println("Check fullcode dummyMerPos: " + dummyMerPos);
            long dummyMer = dummyMerPos>>28;
            long dummyPos = dummyMerPos&mask2;
            
            // Reconstruct (compliment DNA sequence)
//            System.out.println("Check mer sequemce befor compliment : " + dummyMer);
            String binaryMer = Long.toBinaryString(dummyMer);
            int kmer = binaryMer.length()/2;
//            System.out.println("Create compliment at : " + i);
//            System.out.println("Check binaryMer : " + binaryMer);
//            System.out.println("Check dummyPos : " + dummyPos);
            String revBin = new StringBuilder(binaryMer).reverse().toString(); //   reverse Sequence Ex 1011001 to 1001101 
//            System.out.println("Check revBin : " + revBin);
            long revNum = new BigInteger(revBin,2).longValue(); //  Cast binary string to decimal number
            long dummyNewMer = (~revNum)&mask36Bit; //  Create compliment of it
            
//            String strMer = SequenceUtil.decodeMer(dummyMer,kmer);
//            
//            String invMer = SequenceUtil.inverseSequence(strMer);
//            String compMer = SequenceUtil.createComplimentV2(invMer);
           
            //long dummyNewMer = SequenceUtil.encodeMer(compMer, kmer);
            
            //long dummyNewMer = (~dummyMer)&mask36Bit;
//            System.out.println("Check mer sequemce after compliment : " + dummyNewMer);
//            System.out.println("Check Position before reverse : " + dummyPos);
            long dummyNewPos = (this.mersComp.length-1)-dummyPos; /* length-1 because assume array has 10 member ; length is 10 but maximum it index is 9 becaus index start at 0 */
            /* To get new inverse of position value, use this fomular (max index - old index) **Ex. old index is 9 so the inverse of it is (9 - 9) = 0 that's correct!! */
            // Replace to long[] 
//            System.out.println("Check Position after reverse : " + dummyNewPos);
//            System.out.println("Check max index is : " + (this.mers.length-1));
            long dummyNewMerPos = (dummyNewMer<<28)+dummyNewPos;
//            System.out.println("Check fullcode dummyNewMerPost : " + dummyNewMerPos);
            
            if (i%10000000==0){
                System.out.println("Create compliment at : " + i);
            }
            this.mersComp[i] = dummyNewMerPos;
        }
        
        // re-sorted long[]
        Arrays.sort(this.mersComp);
        //return mersComp;
        
    }
    
    public void repeatScan(int iniIdx){
        ArrayList<Integer> listJumpIdx = new ArrayList();
        
        int nextIndex = this.linkIndex[iniIdx];
        long nextMerPos = this.repeatMarkerIndex[nextIndex];
        long nextMer = nextMerPos&this.mask;

        for(int n=this.mainIndex;n<(this.subSequence.length()-numMer)+1;n++){ 
            String sub = this.subSequence.substring(n, n+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
            long compareMer = SequenceUtil.encodeMer(sub, numMer);
            compareMer = compareMer<<28;

            if(nextMer != compareMer){
                break;
            }
        } 
    }
               
               

}
