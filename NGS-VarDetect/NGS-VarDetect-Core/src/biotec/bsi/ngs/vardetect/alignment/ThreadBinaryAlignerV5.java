/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.alignment;

import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.MerRead;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author worawich
 * 
 * Multi-Thread Instruction:
 * 
 * This object will be call in BinaryAlingner on function alignMultithread()
 * 
 * The point that we must separately define new object for multi-thread implement is:
 *  we want to store result or some variable and make it capable to access later 
 *  So we store the result in this object
 *  in order to give the potential of this object to be run in multi-thread 
 *  we Implement Runnable and over ride method run() by put your code tat you want to run it in multi-thread
 *  Next, we have to create some method like start() to create a thread object and start it 
 *  Keep note that when you want to re-run the thread again you have to make sure that the old thread is finish running
 *  by using join() method (the program will stuck at the join() method until all code in run() have done execute)
 *  After that you have to re create the new thread object and start it again.
 *  Don;t forget to access the same object that in case that you have some relative result that you want to store it continuously 
 *  
 * For Version 4 : We change from normal align to align with repeat Marker
 */
public class ThreadBinaryAlignerV5 implements Runnable {
    private Thread t;
    private String threadName;
    private List inputSequence;
    private EncodedSequence encodedRef;
    private long chrNum;
    private int numMer;
    private int threshold;                          // It is a minimum number of count that we accept
    private Map<Long,Long> alignMap;
    private Map<Long,ArrayList<Integer>> alnMerMap;     // Key is align code [strand|alignposition] and value is mer code
    private Map<String,ArrayList<Long>> alnRes;      // Key is ReadName and value is array of long [count|chr|strand|Pos]
    private Map<String,Integer> readLen;            // Key is ReadName and Value is length of each read
    String flag;
            
    
    public ThreadBinaryAlignerV5(String name,List inSeq, EncodedSequence inEncodeRef,long inchr, int inMer , int inThreshold){
        threadName = name;
        inputSequence = inSeq;
        encodedRef = inEncodeRef;
        chrNum = inchr;
        numMer = inMer;
        alnRes = new LinkedHashMap();
        readLen = new LinkedHashMap();
        threshold = inThreshold;
        
        System.out.println("Creating " + threadName);
    }
    
    public void setdata(List inSeq, EncodedSequence inEncodeRef,long inchr, int inMer){    
        inputSequence = inSeq;
        encodedRef = inEncodeRef;
        chrNum = inchr;
        numMer = inMer;
    }
    
    @Override
    public void run(){
        
        System.out.println("Start " + threadName);
        System.out.println("Number of read : " + inputSequence.size());
        /* Alignment algorithm */
        Iterator seqs = inputSequence.iterator();
        while(seqs.hasNext()){                                              // Loop over ShortgunSequence contain in InputSequence 
//            Map<Integer,ArrayList<Integer>> linkIndexCheck = new LinkedHashMap();                       // HashMap contain data that has been use to check for repeat jump
            boolean skipRead = false;
            
            ShortgunSequence seq = (ShortgunSequence)seqs.next();
            Map<Long,ArrayList<Integer>> alnMerMap = new LinkedHashMap();                                         // initialize this hashmap every time when start new loop of Shortgun Read

//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);

            String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence

//                    System.out.print(chr.getName()+" + strand\t");           

            Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
            /* NewPart */
            long oldIniIdx = 0;
            long newIniIdx = 0;
            long iniIndex = 0;
            long recentIdx = 0;
            boolean firstMatchCheck = false;
//            boolean initiateNewReadFlag = true;                                 // this flag is indicate the first time that we consider the read (use to signal inside the Encoded object to renew alnCodeCheckList)
            /************/
            for(int i=0;i<(s.length()-numMer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]
                int index = i;
                String sub = s.substring(i, i+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
                
                if(sub.toUpperCase().equals("AAAAAAAAAAAAAAAAAA")||sub.toUpperCase().equals("TTTTTTTTTTTTTTTTTT")||sub.toUpperCase().equals("GGGGGGGGGGGGGGGGGG")||sub.toUpperCase().equals("CCCCCCCCCCCCCCCCCC")){
                    skipRead = true;
                    break;
                }

//System.out.println("check sub length"+sub.length());
                long m = SequenceUtil.encodeMer(sub, numMer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                
//                if(m == 51904601109L){
//                    System.out.println("");
//                }
                if(m!=-1){                                                          
                    m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose 
//                            long pos = encoded.align(m);
//                    long[] output = encodedRef.align5(m);
                    long[] output = null;
                    /**
                     * Not sure this two line will slow the program or not.
                     */
//                    alnMerMap = (Map<Long,ArrayList<Integer>>)output.get(0);        
//                    alnCodeCheckList = (Map<Long,Long>)output.get(1);
                    /*********/
                            
//                    initiateNewReadFlag = false;
//                    ArrayList<Long> posR = encodedRef.align3(m, s, index, numMer, linkIndexCheck);
//                    if(encodedRef.getRepeatFlag()==false){        
                    if(output == null){   
                        long pos2[] = encodedRef.align5(m);                                // Do alignment with binary search (pos2[] cantain 64 bit long [mer code | position])
    //                            long pos2[] = encoded.fullAlign(m);

                        long pos = -1;
                        if(pos2!=null&&pos2.length>0){
                            pos = pos2[0];
                            pos = pos2.length;
                            //merMap = res.addResult(m, chr.getChrNumber(), pos2);
                        }

                        int idx = (int) (pos-i);
                        if(pos<0){
                          idx = 0;
                        }

                        int totalMer = (seq.getShortgunLength()-numMer)+1;

                        /*************************************************************************************************************/
                        /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/

                        if(pos2 != null){

                            /******** New Part (fixed wrong mer count) Version 3 **********/
                            for(int j=0;j<pos2.length;j++){
                                if(pos2.length>1){
                                    break;
                                }
                                long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition] but already subtract index (offset)


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
                                        iniIndex = alnCodeCheckList.get(alnCode);

                                        indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        merList.remove(0);
                                        merList.add(0,index);
                                        merList.add(1);
                                        alnMerMap.put(indexAlnCode, merList);
                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);

                                        indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                        merList = new ArrayList();
                                        merList.add(0,index);
                                        merList.add(1);
                                        alnMerMap.put(indexAlnCode,merList);
                                    }
                                    /*******************/

                                }else{
                                    iniIndex = index;
                                    alnCodeCheckList.put(alnCode, iniIndex);

                                    long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                    ArrayList<Integer> merList = new ArrayList();
                                    merList.add(0,index);
                                    merList.add(1);
                                    alnMerMap.put(indexAlnCode,merList);

                                }
    //                                    merList = null;
    //                                    System.gc();                                       
                            }

                            /***************************************************************/
                        }
                    }

                    /*-----------------------------------------------------------------------------------------------------------*/
                    /*************************************************************************************************************/
                }
     
            }

            /*************************************************************************************************************/
            /* -------------------------New Implement Part Cons. (Not Stroe in object)---------------------------------------------*/
            
            if(this.alnRes.containsKey(seq.getReadName())&&skipRead==false){                     // Check for existing of ReadName (if exist put result code on existing ArrayList<Long>
                ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                Set keySet = alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                    long count = alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                    long chrIdxStrandAln = (chrNum<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                    long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                    
                    if(count>=threshold){                                            // case check to filter small count peak out (use user specify threshold)
                        countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                    }
                }
                this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
            }else if(this.alnRes.containsKey(seq.getReadName())==false&&skipRead==false){
                ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                Set keySet = alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                    long count = alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                    long chrIdxStrandAln = (chrNum<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                    long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                    
                    if(count>=threshold){                                              // case check to filter small count peak out (use user specify threshold)
                        countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                    }
                }
                this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
            }
            
            this.readLen.put(seq.getReadName(), s.length());
            /* Finish one read clear all data */
//                    this.alnMerMap = null;
//                    System.gc();
            /*-----------------------------------------------------------------------------------------------------------*/
            /*************************************************************************************************************/

//                     System.out.println();
            /* New Implement Part */
            //seq.countAlignmentData(); // Create Alignment count data before change ShortgunSequence

            /*--------------------*/
        }

        /*-------------------- Do compliment alignment -------------------------------*/
        /* Do the same algorithm but use function for compliment */
        Iterator seqsComp = inputSequence.iterator();
//        boolean initiateNewReadFlag = true;                                 // this flag is indicate the first time that we consider the read (use to signal inside the Encoded object to renew alnCodeCheckList)

        while(seqsComp.hasNext()){
//            Map<Integer,ArrayList<Integer>> linkIndexCheck = new LinkedHashMap();                       // HashMap contain data that has been use to check for repeat jump
            boolean skipRead = false;
            ShortgunSequence seq = (ShortgunSequence)seqsComp.next();                                  // get ShortgunSequence from InputSequence
            Map<Long,ArrayList<Integer>> alnMerMap = new LinkedHashMap();

//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);

            String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
            String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
            String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  
//                    System.out.println("******Input Sequence check " + compSeq);
//                    System.out.print(chr.getName()+" - strand\t"); 

            Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
            /* NewPart */
            long oldIniIdx = 0;
            long newIniIdx = 0;
            long iniIndex = 0;
            long recentIdx = 0;
            boolean firstMatchCheck = false;
            /************/
            for(int i=0;i<(compSeq.length()-numMer)+1;i++){                                    // Windowing
                int index = i;                                                              // index at aligncompliment and non compliment is not different. It not effect any thing. we just know strand notation is enough
                String sub = compSeq.substring(i, i+numMer);
                
                if(sub.toUpperCase().equals("AAAAAAAAAAAAAAAAAA")||sub.toUpperCase().equals("TTTTTTTTTTTTTTTTTT")||sub.toUpperCase().equals("GGGGGGGGGGGGGGGGGG")||sub.toUpperCase().equals("CCCCCCCCCCCCCCCCCC")){
                    skipRead = true;
                    break;
                }
                
                //System.out.println("check sub length"+sub.length());
                long m = SequenceUtil.encodeMer(sub, numMer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                if(m!=-1){
                    m = m<<28;
                    
//                    long[] output = encodedRef.align5Compliment(m);
                    long[] output = null;
//                    ArrayList<Long> posR = encodedRef.align3Compliment(m, compSeq, index, numMer, linkIndexCheck);
                    /**
                     * Not sure this two line will slow program or not
                     */
//                    alnMerMap = (Map<Long,ArrayList<Integer>>)output.get(0);        
//                    alnCodeCheckList = (Map<Long,Long>)output.get(1);
                    /***************/

//                    if(encodedRef.getRepeatFlag()==false){
                    if(output == null){ 
                        long pos2[] = encodedRef.align5Compliment(m);                            // Do alignment by alignment function specific for compliment sequence
    //                            long pos2[] = encoded.fullAlign(m);

                        long pos = -1;
                        if(pos2!=null&&pos2.length>0){
                            pos = pos2[0];
                            pos = pos2.length;
                            //merMap = res.addResult(m, chr.getChrNumber(), pos2);
                        }

                        int idx = (int) (pos-i);
                        if(pos<0){
                          idx = 0;
                        }

                        int totalMer = (seq.getShortgunLength()-numMer)+1;

                        /*************************************************************************************************************/
                        /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                                        
                        if(pos2 != null){
                            /******** New Part (fixed wrong mer count) Version 3 **********/
                            for(int j=0;j<pos2.length;j++){
                                if(pos2.length>1){
                                    break;
                                }
                                
                                long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]


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
                                        iniIndex = alnCodeCheckList.get(alnCode);

                                        indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        merList.remove(0);
                                        merList.add(0,index);
                                        merList.add(1);
                                        alnMerMap.put(indexAlnCode, merList);
                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);

                                        indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                        merList = new ArrayList();
                                        merList.add(0,index);
                                        merList.add(1);
                                        alnMerMap.put(indexAlnCode,merList);
                                    }
                                    /*******************/

                                }else{
                                    iniIndex = index;
                                    alnCodeCheckList.put(alnCode, iniIndex);

                                    long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                    ArrayList<Integer> merList = new ArrayList();
                                    merList.add(0,index);
                                    merList.add(1);
                                    alnMerMap.put(indexAlnCode,merList);

                                }
                                       
                            }
                        }
                       
                    }

                    /*-----------------------------------------------------------------------------------------------------------*/
                    /*************************************************************************************************************/
                }
          
            }

            /*************************************************************************************************************/
            /* -------------------------New Implement Part Cons. (Not Stroe in object)---------------------------------------------*/
            
            if(this.alnRes.containsKey(seq.getReadName())&&skipRead==false){                     // Check for existing of ReadName (if exist put result code on existing ArrayList<Long>

                ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                Set keySet = alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                    long count = alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                    long chrIdxStrandAln = (chrNum<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                    long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit

                    if(count>=threshold){                                            // case check to filter small count peak out (use user specify threshold)
                        countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                    }
                }
                this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
            }else if(this.alnRes.containsKey(seq.getReadName())==false && skipRead==false){
                ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                Set keySet = alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                    long count = alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                    long chrIdxStrandAln = (chrNum)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                    long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                    
                    if(count>=threshold){                                                // case check to filter small count peak out (use user specify threshold)
                        countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                    }
                }
                this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
            }
            this.readLen.put(seq.getReadName(), s.length());
            /*-----------------------------------------------------------------------------------------------------------*/
            /*************************************************************************************************************/

            /* New Implement Part */
            //seq.countAlignmentData(); // Create Alignment count data before change ShortgunSequence

            /*--------------------*/
        }
        
        //System.out.println("Thread-"+this.threadName+" : stop");
        this.flag = "run() is Done";
                
    }
    
    public void start(){
        
        t = new Thread (this,threadName);
        t.start();
        this.flag = "run() is running";
        System.out.println("Starting " + threadName + " : " + this.flag);
       
    }
    
    public void join() throws InterruptedException{ 
        t.join();
        System.out.println("Thread-"+this.threadName+" : stop" + " : " + this.flag); 
    }
    
    public Map<String,ArrayList<Long>> getMapResult(){
        return this.alnRes;
    }
    
    public Map<String,Integer> getReadLenList(){
        return this.readLen;
    }
        
}
