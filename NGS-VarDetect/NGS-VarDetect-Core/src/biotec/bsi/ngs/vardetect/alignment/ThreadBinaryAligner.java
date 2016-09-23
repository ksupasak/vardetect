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
 * 
 */
public class ThreadBinaryAligner implements Runnable {
    private Thread t;
    private String threadName;
    private List inputSequence;
    private EncodedSequence encodedRef;
    private long chrNum;
    private int numMer;
    private Map<Long,Long> alignMap;
    private Map<Long,ArrayList<Long>> alnMerMap;     // Key is align code [strand|alignposition] and value is mer code
    private Map<String,ArrayList<Long>> alnRes;      // Key is ReadName and value is array of long [count|chr|strand|Pos]
    String flag;
            
    
    public ThreadBinaryAligner(String name,List inSeq, EncodedSequence inEncodeRef,long inchr, int inMer){
        threadName = name;
        inputSequence = inSeq;
        encodedRef = inEncodeRef;
        chrNum = inchr;
        numMer = inMer;
        alnRes = new LinkedHashMap();
        
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
                    
            ShortgunSequence seq = (ShortgunSequence)seqs.next();
            this.alnMerMap = new HashMap();                                         // initialize this hashmap every time when start new loop of Shortgun Read

//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);

            String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence

//                    System.out.print(chr.getName()+" + strand\t");           


            for(int i=0;i<(s.length()-numMer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]
                int index = i;
                String sub = s.substring(i, i+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
                //System.out.println("check sub length"+sub.length());
                long m = SequenceUtil.encodeMer(sub, numMer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                if(m!=-1){                                                          
                    m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose 
//                            long pos = encoded.align(m);
                    long pos2[] = encodedRef.align2(m);                                // Do alignment with binary search (pos2[] cantain 64 bit long [mer code | position])
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
//                            System.out.println(""+chr.getName()+" "+sub+" "+sub.length()+" : "+m+" pos : "+pos+" idx : "+idx);

//                            System.out.print("\t"+pos);
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("Before Mermap this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            merMap = res.addResult(m, chr.getChrNumber(), pos2); // Still confuse After pss this line all member in pos2 will chang from 28 bit of position to 36 bit of chr:Pos Hoe did it happen !!!
                    // But it work!
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("After mer map this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            res.addResultV2(m, chr.getChrNumber(), pos2, seq.getReadName());
                    //System.out.println("Check seq length" + seq.getShortgunLength());
                    int totalMer = (seq.getShortgunLength()-numMer)+1;

                    /*************************************************************************************************************/
                    /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                    if(pos2 != null){
                        //if(pos2.length == 1){                // (Not work) already check for repeat in same chromosome by checking alignment result must have one match result in this chromosome
                        for(int j=0;j<pos2.length;j++){
                            long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]

                            if(this.alnMerMap.containsKey(alnCode)){
                                ArrayList<Long> merList = this.alnMerMap.get(alnCode);
                                merList.add(m);
                                this.alnMerMap.put(alnCode, merList);    
                            }else{
                                ArrayList<Long> merList = new ArrayList();
                                merList.add(m);
                                this.alnMerMap.put(alnCode,merList);

//                                    merList = null;
//                                    System.gc();
                            }   
                        }    
                    }

                    /*-----------------------------------------------------------------------------------------------------------*/
                    /*************************************************************************************************************/

                }
            //System.out.println(" This mer Map check: "+ (merMap == null));
            //res.createMap(seq.getReadName(), merMap);                    
            }

            /*************************************************************************************************************/
            /* -------------------------New Implement Part Cons. (Not Stroe in object)---------------------------------------------*/

            if(this.alnRes.containsKey(seq.getReadName())){                     // Check for existing of ReadName (if exist put result code on existing ArrayList<Long>
                ArrayList<Long> countChrStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                Set keySet = this.alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long strandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                    long count = this.alnMerMap.get(strandAln).size();          // we can get number of count from number of member in merList
                    long chrStrandAln = (chrNum<<29)+strandAln;     // shift left 29 bit beacause we want to add count number on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<34)+chrStrandAln;          // shift left 34 bit beacause we want to add count number on the front of chrStrandAln which has 34 bit 
                    countChrStrandAlnList.add(countChrStrandAln);
                }
                this.alnRes.put(seq.getReadName(), countChrStrandAlnList);
            }else{
                ArrayList<Long> countChrStrandAlnList = new ArrayList();
                Set keySet = this.alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long strandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                    long count = this.alnMerMap.get(strandAln).size();          // we can get number of count from number of member in merList
                    long chrStrandAln = (chrNum<<29)+strandAln;     // shift left 29 bit beacause we want to add count number on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<34)+chrStrandAln;          // shift left 34 bit beacause we want to add count number on the front of chrStrandAln which has 34 bit 
                    countChrStrandAlnList.add(countChrStrandAln);
                }
                this.alnRes.put(seq.getReadName(), countChrStrandAlnList);
            }


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
        while(seqsComp.hasNext()){
            ShortgunSequence seq = (ShortgunSequence)seqsComp.next();                                  // get ShortgunSequence from InputSequence
            this.alnMerMap = new HashMap();

//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);

            String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
            String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
            String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  
//                    System.out.println("******Input Sequence check " + compSeq);
//                    System.out.print(chr.getName()+" - strand\t"); 


            for(int i=0;i<(compSeq.length()-numMer)+1;i++){                                    // Windowing
                int index = i;                                                              // index at aligncompliment and non compliment is not different. It not effect any thing. we just know strand notation is enough
                String sub = compSeq.substring(i, i+numMer);
                //System.out.println("check sub length"+sub.length());
                long m = SequenceUtil.encodeMer(sub, numMer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                if(m!=-1){
                    m = m<<28;
//                            long pos = encoded.align(m);
                    long pos2[] = encodedRef.align2ComplimentV2(m);                            // Do alignment by alignment function specific for compliment sequence
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
//                            System.out.println(""+chr.getName()+" "+sub+" "+sub.length()+" : "+m+" pos : "+pos+" idx : "+idx);

//                            System.out.print("\t"+pos);
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("Before Mermap this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            merMap = res.addResult(m, chr.getChrNumber(), pos2); // Still confuse After pss this line all member in pos2 will chang from 28 bit of position to 36 bit of chr:Pos Hoe did it happen !!!
                    // But it work!
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("After mer map this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            res.addResultV2(m, chr.getChrNumber(), pos2, seq.getReadName());
                    //System.out.println("Check seq length" + seq.getShortgunLength());
                    int totalMer = (seq.getShortgunLength()-numMer)+1;

                    /*************************************************************************************************************/
                    /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                    if(pos2 != null){
//                                if(pos2.length == 1){           // (Not work) already check for repeat in same chromosome by checking alignment result must have one match result in this chromosome
                        for(int j=0;j<pos2.length;j++){
                            long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]

                            if(this.alnMerMap.containsKey(alnCode)){
                                ArrayList<Long> merList = this.alnMerMap.get(alnCode);
                                merList.add(m);
                                this.alnMerMap.put(alnCode, merList);    
                            }else{
                                ArrayList<Long> merList = new ArrayList();
                                merList.add(m);
                                this.alnMerMap.put(alnCode,merList);

//                                    merList = null;
//                                    System.gc();
                            }
                        }    
                    }

                    /*-----------------------------------------------------------------------------------------------------------*/
                    /*************************************************************************************************************/
                }
            //System.out.println(" This mer Map check: "+ (merMap == null));
            //res.createMap(seq.getReadName(), merMap);                    
            }

            /*************************************************************************************************************/
            /* -------------------------New Implement Part Cons. (Not Stroe in object)---------------------------------------------*/

            if(this.alnRes.containsKey(seq.getReadName())){                     // Check for existing of ReadName (if exist put result code on existing ArrayList<Long>

                ArrayList<Long> countChrStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                Set keySet = this.alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long strandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                    long count = this.alnMerMap.get(strandAln).size();          // we can get number of count from number of member in merList
                    long chrStrandAln = (chrNum<<29)+strandAln;     // shift left 29 bit beacause we want to add count number on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<34)+chrStrandAln;          // shift left 34 bit beacause we want to add count number on the front of chrStrandAln which has 34 bit

                    countChrStrandAlnList.add(countChrStrandAln);
                }
                this.alnRes.put(seq.getReadName(), countChrStrandAlnList);
            }else{
                ArrayList<Long> countChrStrandAlnList = new ArrayList();
                Set keySet = this.alnMerMap.keySet();
                Iterator keyIter =keySet.iterator();
                while(keyIter.hasNext()){
                    long strandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                    long count = this.alnMerMap.get(strandAln).size();          // we can get number of count from number of member in merList
                    long chrStrandAln = (chrNum<<29)+strandAln;     // shift left 29 bit beacause we want to add count number on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<34)+chrStrandAln;          // shift left 34 bit beacause we want to add count number on the front of chrStrandAln which has 34 bit 
                    countChrStrandAlnList.add(countChrStrandAln);
                }
                this.alnRes.put(seq.getReadName(), countChrStrandAlnList);
            }


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
        
}
