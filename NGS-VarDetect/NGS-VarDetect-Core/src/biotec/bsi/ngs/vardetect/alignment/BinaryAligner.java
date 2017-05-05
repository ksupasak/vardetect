/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.alignment;

import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResult;
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.MerRead;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author soup
 */

public class BinaryAligner extends Thread implements Aligner {

    ReferenceSequence ref;
    int mer = 18;
    int numberOfThread = 1;
    ArrayList<ThreadBinaryAligner> threadList;
    ArrayList<ThreadBinaryAlignerV3> threadListV3;
    ArrayList<ThreadBinaryAlignerV4> threadListV4;
    ArrayList<ThreadBinaryAlignerV5> threadListV5;
    ArrayList<ThreadBinaryAlignerLongRead> threadListLongRead;
    Map<Long,ArrayList<Long>> alnMerMap;     // Key is align code [strand|alignposition] and value is mer code
    Map<String,ArrayList<Long>> alnRes;      // Key is ReadName and value is array of long [count|chr|iniIdx|strand|Pos]
    Map<String,ArrayList<Long>> alnRes1;      // Key is ReadName and value is array of long [iniIndex|strand|Pos]
    Map<String,ArrayList<Long>> alnRes2;      // Key is ReadName and value is array of long [chr|count]
    Map<String,Integer> readLen;      // Key is ReadName and value is Integer of read length    
    Map<Long,Byte> varType;
    Map<String,ArrayList<Byte>> varBox;
    
    public void setReferenceSequence(ReferenceSequence ref){
        this.ref = ref;
        
    }
       
    public AlignmentResultRead align(ReferenceSequence ref, InputSequence input, int numMer, int threshold) {
        
        this.setReferenceSequence(ref);
        return align(input,numMer,threshold);
        
    }
    
    public AlignmentResultRead align(InputSequence input, int numMer, int threshold){
        
        
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;
        
        
        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            try {
                
                
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = SequenceUtil.createAllReference(chr, this.mer);   // Create or import all reference [chromosome reference, repeat index, repeat Marker]
                
                while(seqs.hasMoreElements()){                                              // Loop over ShortgunSequence contain in InputSequence 
                    ShortgunSequence seq = seqs.nextElement();
                    Map<Long,long[]> merMap = new HashMap();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence
                    
//                    System.out.print(chr.getName()+" + strand\t");           
                    
                    
                    for(int i=0;i<(s.length()-mer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]
                        int index = i;
                        String sub = s.substring(i, i+mer);                                 // cut String sequence into sub string sequence (mer length long) 
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){                                                          
                            m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose 
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2(m);                                // Do alignment with binary search (pos2[] cantain 29 bit long [strand notation| position])
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Stroe in object)---------------------------------------------*/
                            //if(pos2 != null){
                            if (seq.getMerReadSize() < totalMer){                                   // Check selected ShortgunSequence has all possible caontain in it or not 
                                MerRead merRead = new MerRead();                                    // Create merRead object                                
                                merRead.addMatchResultStrand(m, pos2, index, chr.getChrNumber(),this.mer);   // add mer code(36 bit), align result(64 bit [mer|pos]), index and chromosome number into MerRead
                                seq.addMerRead(merRead);                                            // add object MerRead into ShortgunSequence
//                                System.out.println("(First Time) Size Mer Read check: " + seq.getMerReadSize());
                            }else{
//                                System.out.println("(Other time) Size Mer Read check: " + seq.getMerReadSize());
                                if (index<seq.getMerReadSize()){
//                                    System.out.println("Check Index: " + index);
                                    //seq.getMerRead().contains()
                                    MerRead merRead = seq.getMerRead().get(index);                  // call back existing merRead to update                                   
                                    merRead.addMatchResultStrand(m, pos2, index, chr.getChrNumber(),this.mer); // add mer code(36 bit), align result(64 bit [mer|pos]), index and chromosome number into MerRead 
                                    seq.addMerReadByIndex(index,merRead);                           // add MerRead into ShortgunSequence by index
                                    
//                                    System.out.println("(Other time) Should be constant: Size Mer Read check: " + seq.getMerReadSize());
                                }                                   
                            }
                            //}
                            /*-----------------------------------------------------------------------------------------------------------*/
                            /*************************************************************************************************************/
                            
                        }
                    //System.out.println(" This mer Map check: "+ (merMap == null));
                    //res.createMap(seq.getReadName(), merMap);                    
                    }
//                     System.out.println();
                    /* New Implement Part */
                    //seq.countAlignmentData(); // Create Alignment count data before change ShortgunSequence
                    
                    /*--------------------*/
                }
                
                /*-------------------- Do compliment alignment -------------------------------*/
                /* Do the same algorithm but use function for compliment */
                Enumeration<ShortgunSequence> seqsComp = input.getInputSequence().elements();
                while(seqsComp.hasMoreElements()){
                    ShortgunSequence seq = seqsComp.nextElement();                                  // get ShortgunSequence from InputSequence
                    Map<Long,long[]> merMap = new HashMap();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
                    String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  
//                    System.out.println("******Input Sequence check " + compSeq);
//                    System.out.print(chr.getName()+" - strand\t"); 
                    
                    
                    for(int i=0;i<(compSeq.length()-mer)+1;i++){                                    // Windowing
                        int index = i;
                        String sub = compSeq.substring(i, i+mer);
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){
                            m = m<<28;
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2ComplimentV2(m);                            // Do alignment by alignment function specific for compliment sequence
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Stroe in object)---------------------------------------------*/
                            if (seq.getMerReadSize() < totalMer){
                                MerRead merRead = new MerRead();
                                
                                merRead.addMatchResultStrandCompliment(m, pos2, index, chr.getChrNumber(),this.mer);             // add infomation into MerRead
                                seq.addMerRead(merRead);
                                System.out.println("*********** This word should not apear ***************");
                            }else{
//                                System.out.println("(Other time) Size Mer Read check: " + seq.getMerReadSize());
                                if (index<seq.getMerReadSize()){
                                    int compIndex = (seq.getMerReadSize()-1) - index;       // not use (does not effect any thing)  
//                                    System.out.println("Check Index: " + index);
                                    MerRead merRead = seq.getMerRead().get(index);
                                    merRead.addMatchResultStrandCompliment(m, pos2, index, chr.getChrNumber(),this.mer);
                                    seq.addMerReadByIndex(index,merRead);
//                                    MerRead merRead = seq.getMerRead().get(compIndex); // call back existing merRead to update                                    
//                                    merRead.addMatchResultStrand(pos2, chr.getChrNumber());
//                                    seq.addMerReadByIndex(compIndex,merRead);
//                                    System.out.println("(Other time) Should be constant: Size Mer Read check: " + seq.getMerReadSize());
                                }                                   
                            }
                            /*-----------------------------------------------------------------------------------------------------------*/
                            /*************************************************************************************************************/
                            
                        }
                    //System.out.println(" This mer Map check: "+ (merMap == null));
                    //res.createMap(seq.getReadName(), merMap);                    
                    }
//                     System.out.println();
                    /* New Implement Part */
                    //seq.countAlignmentData(); // Create Alignment count data before change ShortgunSequence
                    
                    /*--------------------*/
                }
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();
                
                
                
            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        /* After Alignment we will loop over ShortgunSequence and do the Alignmentcount */
        //AlignmentResultRead alinResult = new AlignmentResultRead();  
        Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
        while(seqs.hasMoreElements()){
                    ShortgunSequence seq = seqs.nextElement();
                    seq.countAlignmentData();                                       // count alignment data (this function must be call after information is ready)
                    alinResult.addResult(seq);                                      // add ShortgunSequence into AlignmentResultRead  
        }
        
        //return res;
        return alinResult;
    }
    
    public AlignmentResultRead alignMultithread(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException {
        
        this.setReferenceSequence(ref);
        return alignMultithread(input,numThread,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignMultithread(InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException{
        
        /**
        * This method will create object that implement multi-thread capability in it
        * This function will split input sequence into number of portion relate to number of user specify thread number 
        * Then it will loop to create a set of thread and store in threadList
        * Each portion of split input will be pass into each thread implement object and do there job
        * At last we will wait for each thread finish there job and merge the result together
        * 
        */
        
        threadList = new ArrayList();
        this.numberOfThread = numThread;
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;

        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            try {
                
                int count = 0;
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = SequenceUtil.createAllReference(chr, this.mer);   // Create or import all reference [chromosome reference, repeat index, repeat Marker]
                long chrnumber = chr.getChrNumber();
                /*********/
                int inputSize = input.getInputSequence().size();
                double dummyNum = (double)inputSize/this.numberOfThread;
                int numPerPartition = (int)Math.ceil(dummyNum);
                
                /* Create thread object follow by specific numThread */
                /*  Separate case In order to use existing thread object or create new for first time*/
                
                if(threadList.isEmpty()){
                 
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){                        
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        String threadName = "Thread_"+count; 
                        ThreadBinaryAligner newThread = new ThreadBinaryAligner(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer);
                        newThread.start();
                        threadList.add(newThread);
                        count++;
                    }
                }else{
                    int dummynum = 0;
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        ThreadBinaryAligner dummythread = threadList.get(dummynum);
                        dummythread.setdata(splitInputSequence, encoded, chrnumber, mer);
                        dummythread.start();
                        dummynum++;
                    }
                }
                              
//                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
//                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
//                        String threadName = "Thread_"+count; 
//                        ThreadBinaryAligner newThread = new ThreadBinaryAligner(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer);
//                        newThread.start();
//
//                        threadList.add(newThread);
//                        count++;
//                   
                              
                System.out.println("Number of thread check : " + threadList.size() + " Must equal to " + numThread);
                for(int i=0;i<threadList.size();i++){
                    /* Wait for specific thread to finish execute (finish method run()) */ 
                    /* After run() method is finish thread stop execute but the thread object that we overwrite run() in it is still exist */
                    threadList.get(i).join();                    
                }
                
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();

            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        
        for(int i=0;i<threadList.size();i++){
            /* Loop to get align result map from each thread object */
            /* Put all align map result into one Map result */
            if(i==0){
                this.alnRes = threadList.get(i).getMapResult();
            }else{
                this.alnRes.putAll(threadList.get(i).getMapResult());
            }
        }
        
        
        alinResult.addMapResult(this.alnRes);
 
        return alinResult;
    }
    
    public AlignmentResultRead alignMultithreadV3(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException {
        /**
         * This version 3 function give the v3 data structure of result
         */
        this.setReferenceSequence(ref);
        return alignMultithreadV3(input,numThread,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignMultithreadV3(InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException{
        
        /**
        * This method will create object that implement multi-thread capability in it
        * This function will split input sequence into number of portion relate to number of user specify thread number 
        * Then it will loop to create a set of thread and store in threadList
        * Each portion of split input will be pass into each thread implement object and do there job
        * At last we will wait for each thread finish there job and merge the result together
        * 
        * This version 3 function give the v3 data structure of result
        */
        
        threadListV3 = new ArrayList();
        this.numberOfThread = numThread;
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;

        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            long startTime = System.currentTimeMillis();
            try {
                
                int count = 0;
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = SequenceUtil.createAllReference(chr, this.mer);   // Create or import all reference [chromosome reference, repeat index, repeat Marker]
                long chrnumber = chr.getChrNumber();
                /*********/
                int inputSize = input.getInputSequence().size();
                double dummyNum = (double)inputSize/this.numberOfThread;
                int numPerPartition = (int)Math.ceil(dummyNum);
                
                /* Create thread object follow by specific numThread */
                /*  Separate case In order to use existing thread object or create new for first time*/
                
                if(threadListV3.isEmpty()){
                 
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){                        
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        String threadName = "Thread_"+count; 
                        ThreadBinaryAlignerV3 newThread = new ThreadBinaryAlignerV3(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer,threshold);
                        newThread.start();                        
                        threadListV3.add(newThread);
                        count++;
                    }
                }else{
                    int dummynum = 0;
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        ThreadBinaryAlignerV3 dummythread = threadListV3.get(dummynum);
                        dummythread.setdata(splitInputSequence, encoded, chrnumber, mer);
                        dummythread.start();
                        dummynum++;
                    }
                }
                              
//                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
//                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
//                        String threadName = "Thread_"+count; 
//                        ThreadBinaryAligner newThread = new ThreadBinaryAligner(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer);
//                        newThread.start();
//
//                        threadList.add(newThread);
//                        count++;
//                   
                              
                System.out.println("Number of thread check : " + threadListV3.size() + " Must equal to " + numThread);
                for(int i=0;i<threadListV3.size();i++){
                    /* Wait for specific thread to finish execute (finish method run()) */ 
                    /* After run() method is finish thread stop execute but the thread object that we overwrite run() in it is still exist */
                    threadListV3.get(i).join();                    
                }
                
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();

            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            long endTime = System.currentTimeMillis();
            double totalTime = (endTime - startTime)/1000;
            System.out.println("Done time use: "+ totalTime +" second");
            System.out.println();
        }
        
        for(int i=0;i<threadListV3.size();i++){
            /* Loop to get align result map from each thread object */
            /* Put all align map result into one Map result */
            if(i==0){
                this.alnRes = threadListV3.get(i).getMapResult();
            }else{
                this.alnRes.putAll(threadListV3.get(i).getMapResult());
            }
        }
        
        
        alinResult.addMapResult(this.alnRes);
 
        return alinResult;
    }
    
    public AlignmentResultRead alignMultithreadV4(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException {
        /**
         * This version 3 function give the v3 data structure of result
         */
        this.setReferenceSequence(ref);
        return alignMultithreadV4(input,numThread,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignMultithreadV4(InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException{
        
        /**
        * This method will create object that implement multi-thread capability in it
        * This function will split input sequence into number of portion relate to number of user specify thread number 
        * Then it will loop to create a set of thread and store in threadList
        * Each portion of split input will be pass into each thread implement object and do there job
        * At last we will wait for each thread finish there job and merge the result together
        * 
        * This version 3 function give the v3 data structure of result
        */
        
        threadListV4 = new ArrayList();
        this.numberOfThread = numThread;
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;

        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            long startTime = System.currentTimeMillis();
            try {
                
                int count = 0;
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");
                EncodedSequence encoded = SequenceUtil.createAllReferenceV2(chr, this.mer, 'a');   // Create or import all reference [chromosome reference, repeat index, repeat Marker]                               
                long chrnumber = chr.getChrNumber();
                /*********/
                int inputSize = input.getInputSequence().size();
                double dummyNum = (double)inputSize/this.numberOfThread;
                int numPerPartition = (int)Math.ceil(dummyNum);
                
                /* Create thread object follow by specific numThread */
                /*  Separate case In order to use existing thread object or create new for first time*/
                
                if(threadListV4.isEmpty()){
                 
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){                        
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        String threadName = "Thread_"+count; 
                        ThreadBinaryAlignerV4 newThread = new ThreadBinaryAlignerV4(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer,threshold);
                        newThread.start();                        
                        threadListV4.add(newThread);
                        count++;
                    }
                }else{
                    int dummynum = 0;
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        ThreadBinaryAlignerV4 dummythread = threadListV4.get(dummynum);
                        dummythread.setdata(splitInputSequence, encoded, chrnumber, mer);
                        dummythread.start();
                        dummynum++;
                    }
                }
                              
//                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
//                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
//                        String threadName = "Thread_"+count; 
//                        ThreadBinaryAligner newThread = new ThreadBinaryAligner(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer);
//                        newThread.start();
//
//                        threadList.add(newThread);
//                        count++;
//                   
                              
                System.out.println("Number of thread check : " + threadListV4.size() + " Must equal to " + numThread);
                for(int i=0;i<threadListV4.size();i++){
                    /* Wait for specific thread to finish execute (finish method run()) */ 
                    /* After run() method is finish thread stop execute but the thread object that we overwrite run() in it is still exist */
                    threadListV4.get(i).join();                    
                }
                
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();

            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            long endTime = System.currentTimeMillis();
            double totalTime = (endTime - startTime)/1000;
            System.out.println("Done time use: "+ totalTime +" second");
            System.out.println();
        }
        
        for(int i=0;i<threadListV4.size();i++){
            /* Loop to get align result map from each thread object */
            /* Put all align map result into one Map result */
            if(i==0){
                this.alnRes = threadListV4.get(i).getMapResult();
                this.readLen = threadListV4.get(i).getReadLenList();
            }else{
                this.alnRes.putAll(threadListV4.get(i).getMapResult());
                this.readLen.putAll(threadListV4.get(i).getReadLenList());
            }
        }
        
        
        alinResult.addMapResult(this.alnRes);
        alinResult.addReadLenList(this.readLen);
 
        return alinResult;
    }
    
    public AlignmentResultRead alignMultithreadV5(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException {
        /**
         * This version 3 function give the v3 data structure of result
         * compatible with short read
         * 
         * already implement read length transfer  
         */
        this.setReferenceSequence(ref);
        return alignMultithreadV5(input,numThread,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignMultithreadV5(InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException{
        
        /**
        * This method will create object that implement multi-thread capability in it
        * This function will split input sequence into number of portion relate to number of user specify thread number 
        * Then it will loop to create a set of thread and store in threadList
        * Each portion of split input will be pass into each thread implement object and do there job
        * At last we will wait for each thread finish there job and merge the result together
        * 
        * This version 3 function give the v3 data structure of result
        * Conpatible with short Read
        * Limit read length 255 bp (8bit)
        * Limit number of chromosome 31 chromosome (5bit)
        * 
        * this version has build for cut repeat protocol (ignore all repeat consider only unique in each chromosome)
        * already implement read length transfer
        */
        
        threadListV5 = new ArrayList();
        this.numberOfThread = numThread;
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;

        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            long startTime = System.currentTimeMillis();
            try {
                
                int count = 0;
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");
                EncodedSequence encoded = SequenceUtil.createAllReferenceV2(chr, this.mer, 'r');   // Create or import all reference [chromosome reference, repeat index, repeat Marker]                               
                long chrnumber = chr.getChrNumber();
                /*********/
                int inputSize = input.getInputSequence().size();
                double dummyNum = (double)inputSize/this.numberOfThread;
                int numPerPartition = (int)Math.ceil(dummyNum);
                
                /* Create thread object follow by specific numThread */
                /*  Separate case In order to use existing thread object or create new for first time*/
                
                if(threadListV5.isEmpty()){
                 
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){                        
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        String threadName = "Thread_"+count; 
                        ThreadBinaryAlignerV5 newThread = new ThreadBinaryAlignerV5(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer,threshold);
                        newThread.start();                        
                        threadListV5.add(newThread);
                        count++;
                    }
                }else{
                    int dummynum = 0;
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        ThreadBinaryAlignerV5 dummythread = threadListV5.get(dummynum);
                        dummythread.setdata(splitInputSequence, encoded, chrnumber, mer);
                        dummythread.start();
                        dummynum++;
                    }
                }
                              
//                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
//                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
//                        String threadName = "Thread_"+count; 
//                        ThreadBinaryAligner newThread = new ThreadBinaryAligner(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer);
//                        newThread.start();
//
//                        threadList.add(newThread);
//                        count++;
//                   
                              
                System.out.println("Number of thread check : " + threadListV5.size() + " Must equal to " + numThread);
                for(int i=0;i<threadListV5.size();i++){
                    /* Wait for specific thread to finish execute (finish method run()) */ 
                    /* After run() method is finish thread stop execute but the thread object that we overwrite run() in it is still exist */
                    threadListV5.get(i).join();                    
                }
                
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();

            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            long endTime = System.currentTimeMillis();
            double totalTime = (endTime - startTime)/1000;
            System.out.println("Done time use: "+ totalTime +" second");
            System.out.println();
        }
        
        for(int i=0;i<threadListV5.size();i++){
            /* Loop to get align result map from each thread object */
            /* Put all align map result into one Map result */
            if(i==0){
                this.alnRes = threadListV5.get(i).getMapResult();
                this.readLen = threadListV5.get(i).getReadLenList();
            }else{
                this.alnRes.putAll(threadListV5.get(i).getMapResult());
                this.readLen.putAll(threadListV5.get(i).getReadLenList());
            }
        }
        
        
        alinResult.addMapResult(this.alnRes);
        alinResult.addReadLenList(this.readLen);
 
        return alinResult;
    }
    
    public AlignmentResultRead alignMultithreadLongRead(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException {
        /**
         * This version 3 function give the v3 data structure of result
         * Compatible with long read 
         */
        this.setReferenceSequence(ref);
        return alignMultithreadLongRead(input,numThread,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignMultithreadLongRead(InputSequence input, int numThread, int numMer, int threshold) throws InterruptedException{
        
        /**
        * This method will create object that implement multi-thread capability in it
        * This function will split input sequence into number of portion relate to number of user specify thread number 
        * Then it will loop to create a set of thread and store in threadList
        * Each portion of split input will be pass into each thread implement object and do there job
        * At last we will wait for each thread finish there job and merge the result together
        * 
        * This version 3 function give the v3 data structure of result
        * 
        * this version has build for cut repeat protocol (ignore all repeat consider only unique in each chromosome)
        * Compatible with longRead
        * Limit max chromosome at 31 chromosome (5bit)
        */
        
        threadListLongRead = new ArrayList();
        this.numberOfThread = numThread;
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;

        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            long startTime = System.currentTimeMillis();
            try {
                
                int count = 0;
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");
                EncodedSequence encoded = SequenceUtil.createAllReferenceV2(chr, this.mer, 'r');   // Create or import all reference [chromosome reference, repeat index, repeat Marker]                               
                long chrnumber = chr.getChrNumber();
                /*********/
                int inputSize = input.getInputSequence().size();
                double dummyNum = (double)inputSize/this.numberOfThread;
                int numPerPartition = (int)Math.ceil(dummyNum);
                
                /* Create thread object follow by specific numThread */
                /*  Separate case In order to use existing thread object or create new for first time*/
                
                if(threadListLongRead.isEmpty()){
                 
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){                        
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        String threadName = "Thread_"+count; 
                        ThreadBinaryAlignerLongRead newThread = new ThreadBinaryAlignerLongRead(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer,threshold);
                        newThread.start();                        
                        threadListLongRead.add(newThread);
                        count++;
                    }
                }else{
                    int dummynum = 0;
                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
                        ThreadBinaryAlignerLongRead dummythread = threadListLongRead.get(dummynum);
                        dummythread.setdata(splitInputSequence, encoded, chrnumber, mer);
                        dummythread.start();
                        dummynum++;
                    }
                }
                              
//                    for(int i = 0 ; i < inputSize ; i+= numPerPartition){
//                        List splitInputSequence = input.getInputSequence().subList(i, Math.min(inputSize, i+numPerPartition));
//                        String threadName = "Thread_"+count; 
//                        ThreadBinaryAligner newThread = new ThreadBinaryAligner(threadName,splitInputSequence,encoded,chr.getChrNumber(),mer);
//                        newThread.start();
//
//                        threadList.add(newThread);
//                        count++;
//                   
                              
                System.out.println("Number of thread check : " + threadListLongRead.size() + " Must equal to " + numThread);
                for(int i=0;i<threadListLongRead.size();i++){
                    /* Wait for specific thread to finish execute (finish method run()) */ 
                    /* After run() method is finish thread stop execute but the thread object that we overwrite run() in it is still exist */
                    threadListLongRead.get(i).join();                    
                }
                
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();

            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            long endTime = System.currentTimeMillis();
            double totalTime = (endTime - startTime)/1000;
            System.out.println("Done time use: "+ totalTime +" second");
            System.out.println();
        }
        
        for(int i=0;i<threadListLongRead.size();i++){
            /* Loop to get align result map from each thread object */
            /* Put all align map result into one Map result */
            if(i==0){
                this.alnRes1 = threadListLongRead.get(i).getMapResult1();
                this.alnRes2 = threadListLongRead.get(i).getMapResult2();
                this.readLen = threadListLongRead.get(i).getReadLenList();
            }else{
                this.alnRes1.putAll(threadListLongRead.get(i).getMapResult1());
                this.alnRes2.putAll(threadListLongRead.get(i).getMapResult2());
                this.readLen.putAll(threadListLongRead.get(i).getReadLenList());
            }
        }       
        
        
        alinResult.addMapResult1(this.alnRes1);
        alinResult.addMapResult2(this.alnRes2);
        alinResult.addReadLenList(this.readLen);
        return alinResult;
    }
    
    public AlignmentResultRead alignV2(ReferenceSequence ref, InputSequence input, int numMer, int threshold) {
        
        this.setReferenceSequence(ref);
        return alignV2(input,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignV2(InputSequence input, int numMer, int threshold){
        
        alnRes = new LinkedHashMap();                                     // initialize this hashmap one time per alinment job
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;
        
        
        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            try {
                
                
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = SequenceUtil.createAllReference(chr, this.mer);   // Create or import all reference [chromosome reference, repeat index, repeat Marker]
                while(seqs.hasMoreElements()){                                              // Loop over ShortgunSequence contain in InputSequence 
                    ShortgunSequence seq = seqs.nextElement();
                    this.alnMerMap = new HashMap();                                         // initialize this hashmap every time when start new loop of Shortgun Read
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence
                    
//                    System.out.print(chr.getName()+" + strand\t");           
                    
                    
                    for(int i=0;i<(s.length()-mer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]
                        int index = i;
                        String sub = s.substring(i, i+mer);                                 // cut String sequence into sub string sequence (mer length long) 
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){                                                          
                            m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose 
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2(m);                                // Do alignment with binary search (pos2[] cantain 64 bit long [mer code | position])
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
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
                            long chrStrandAln = (chr.getChrNumber()<<29)+strandAln;     // shift left 29 bit beacause we want to add count number on the front of strandAln which has 29 bit
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
                            long chrStrandAln = (chr.getChrNumber()<<29)+strandAln;     // shift left 29 bit beacause we want to add count number on the front of strandAln which has 29 bit
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
                Enumeration<ShortgunSequence> seqsComp = input.getInputSequence().elements();
                while(seqsComp.hasMoreElements()){
                    ShortgunSequence seq = seqsComp.nextElement();                                  // get ShortgunSequence from InputSequence
                    this.alnMerMap = new HashMap();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
                    String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  
//                    System.out.println("******Input Sequence check " + compSeq);
//                    System.out.print(chr.getName()+" - strand\t"); 
                    
                    
                    for(int i=0;i<(compSeq.length()-mer)+1;i++){                                    // Windowing
                        int index = i;                                                              // index at aligncompliment and non compliment is not different. It not effect any thing. we just know strand notation is enough
                        String sub = compSeq.substring(i, i+mer);
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){
                            m = m<<28;
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2ComplimentV2(m);                            // Do alignment by alignment function specific for compliment sequence
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
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
                            long chrStrandAln = (chr.getChrNumber()<<29)+strandAln;     // shift left 29 bit beacause we want to add chr number on the front of strandAln which has 29 bit
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
                            long chrStrandAln = (chr.getChrNumber()<<29)+strandAln;     // shift left 29 bit beacause we want to add chr number on the front of strandAln which has 29 bit
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
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();
                
                
                
            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        alinResult.addMapResult(this.alnRes);
        /* After Alignment we will loop over ShortgunSequence and do the Alignmentcount */
        //AlignmentResultRead alinResult = new AlignmentResultRead();  
//        Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
//        while(seqs.hasMoreElements()){
//                    ShortgunSequence seq = seqs.nextElement();
//                    seq.countAlignmentData();                                       // count alignment data (this function must be call after information is ready)
//                    alinResult.addResult(seq);                                      // add ShortgunSequence into AlignmentResultRead  
//        }
        
        //return res;
        return alinResult;
    }
    
    public AlignmentResultRead localAlign(EncodedSequence ref, InputSequence input, int kmer, long numberOflocalRef){
        
        mer = kmer;
        alnRes = new LinkedHashMap();                                     // initialize this hashmap one time per alinment job
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        
        Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
        EncodedSequence encoded = ref;                                             // encoded selected chromosome (just for sure it is encode)
        while(seqs.hasMoreElements()){                                              // Loop over ShortgunSequence contain in InputSequence 
            ShortgunSequence seq = seqs.nextElement();                              // seq is input Read
            this.alnMerMap = new HashMap();                                         // initialize this hashmap every time when start new loop of Shortgun Read

//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);

            String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence

//                    System.out.print(chr.getName()+" + strand\t");           


            for(int i=0;i<(s.length()-mer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]
                int index = i;
                String sub = s.substring(i, i+mer);                                 // cut String sequence into sub string sequence (mer length long) 
                //System.out.println("check sub length"+sub.length());
                long m = SequenceUtil.encodeMer(sub, mer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                if(m!=-1){                                                          
                    m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose 
//                            long pos = encoded.align(m);
                    long pos2[] = encoded.align2(m);                                // Do alignment with binary search (pos2[] cantain 64 bit long [mer code | position])
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
            
                    int totalMer = (seq.getShortgunLength()-mer)+1;

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
                    long chrStrandAln = (numberOflocalRef<<29)+strandAln;     // shift left 29 bit beacause we want to add numberOflocalRef on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<56)+chrStrandAln;          // shift left 56 bit beacause we want to add count number on the front of chrStrandAln which has 56 from (29+27) bit (we reserve bit for numberOflocalRef for 27 bit) 
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
                    long chrStrandAln = (numberOflocalRef<<29)+strandAln;     // shift left 29 bit beacause we want to add numberOflocalRef on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<56)+chrStrandAln;          // shift left 56 bit beacause we want to add count number on the front of chrStrandAln which has 56 from (29+27) bit (we reserve bit for numberOflocalRef for 27 bit)
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
        Enumeration<ShortgunSequence> seqsComp = input.getInputSequence().elements();
        while(seqsComp.hasMoreElements()){
            ShortgunSequence seq = seqsComp.nextElement();                                  // get ShortgunSequence from InputSequence
            this.alnMerMap = new HashMap();

//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);

            String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
            String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
            String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  
//                    System.out.println("******Input Sequence check " + compSeq);
//                    System.out.print(chr.getName()+" - strand\t"); 


            for(int i=0;i<(compSeq.length()-mer)+1;i++){                                    // Windowing
                int index = i;                                                              // index at aligncompliment and non compliment is not different. It not effect any thing. we just know strand notation is enough
                String sub = compSeq.substring(i, i+mer);
                //System.out.println("check sub length"+sub.length());
                long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                if(m!=-1){
                    m = m<<28;
//                            long pos = encoded.align(m);
                    long pos2[] = encoded.align2ComplimentV2(m);                            // Do alignment by alignment function specific for compliment sequence
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
                    
                    int totalMer = (seq.getShortgunLength()-mer)+1;

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
                    long chrStrandAln = (numberOflocalRef<<29)+strandAln;     // shift left 29 bit beacause we want to add numberOflocalRef on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<56)+chrStrandAln;          // shift left 56 bit beacause we want to add count number on the front of chrStrandAln which has 56 from (29+27) bit (we reserve bit for numberOflocalRef for 27 bit)

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
                    long chrStrandAln = (numberOflocalRef<<29)+strandAln;     // shift left 29 bit beacause we want to add numberOflocalRef on the front of strandAln which has 29 bit
                    long countChrStrandAln = (count<<56)+chrStrandAln;          // shift left 56 bit beacause we want to add count number on the front of chrStrandAln which has 56 from (29+27) bit (we reserve bit for numberOflocalRef for 27 bit)
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
        
        alinResult.addMapResult(this.alnRes);       // alnRes is Map of each input(key) and align resut(value)
        
        return alinResult;
    }

    public AlignmentResultRead alignV3(ReferenceSequence ref, InputSequence input, int numMer, int threshold) {
        
        this.setReferenceSequence(ref);
        return alignV3(input,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignV3(InputSequence input, int numMer, int threshold){
        /**
         * The different between alignV2 and V3 is we add initial match of read index (iniIdx) 
         * So all sorting or other old version function may not suitable for this new data structure
         * Please use the version of it;
         */
        alnRes = new LinkedHashMap();                                     // initialize this hashmap one time per alinment job
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;
        
        
        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            try {
                
                
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = SequenceUtil.createAllReference(chr, this.mer);   // Create or import all reference [chromosome reference, repeat index, repeat Marker]               
                while(seqs.hasMoreElements()){                                              // Loop over ShortgunSequence contain in InputSequence 
                    ShortgunSequence seq = seqs.nextElement();
                    this.alnMerMap = new LinkedHashMap();                                         // initialize this hashmap every time when start new loop of Shortgun Read
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence
                    
//                    System.out.print(chr.getName()+" + strand\t");           
                    
                    Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
                    /******** New Part (fixed wrong mer count) **********/
                    long oldIniIdx = 0;
                    long newIniIdx = 0;
                    long iniIndex = 0;
                    long recentIdx = 0;
                    boolean firstMatchCheck = false;
                    /****************/
                    for(int i=0;i<(s.length()-mer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]    
                        int index = i;
                        String sub = s.substring(i, i+mer);                                 // cut String sequence into sub string sequence (mer length long) 
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){                                                          
                            m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose 
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2(m);                                // Do alignment with binary search (pos2[] cantain 64 bit long [mer code | position])
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                            /*------------------------------- add iniIndex in front of 29 bit [strand|position] -----------------------------*/ 
                            if(pos2 != null){
                                //if(pos2.length == 1){                // (Not work) already check for repeat in same chromosome by checking alignment result must have one match result in this chromosome
                                
                                
                                /* Old Part *//*
                                for(int j=0;j<pos2.length;j++){
                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
                                    
                                    
                                    if(alnCodeCheckList.containsKey(alnCode)){
                                        long iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode, merList);
                                        
                                    }else{
                                        long iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode,merList);
                                        
                                    }
    //                                    merList = null;
    //                                    System.gc();                                       
                                }
                                */
                                
                                /******** New Part (fixed wrong mer count) Version 3 **********/
                                for(int j=0;j<pos2.length;j++){
                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
                                    
                                    
                                    if(alnCodeCheckList.containsKey(alnCode)){
 
                                        iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
                                        
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
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode, merList);
                                        }else{
                                            /**
                                             * it's not continue. So, iniIndex has change to present index                                                                                  
                                             */
                                            
                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode,merList);
                                        }
                                        /*******************/
                                        
                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(0,(long)index);
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode,merList);
                                        
                                    }
    //                                    merList = null;
    //                                    System.gc();                                       
                                }
                                
                                /***************************************************************/
                                
                                
                                
                                /******** New Part (fixed wrong mer count) Version 2 (not work)**********/
                                
                                /* Specify iniIdx check from continueously of index */
//                                newIniIdx = index;
//                                if(firstMatchCheck==false){
//                                    // it's first time
//                                    iniIndex = index;
//                                    oldIniIdx = index;
////                                    recentIdx = index;
//                                    firstMatchCheck = true;
//                                }else{
//                                    if(newIniIdx-recentIdx==1){
//                                        iniIndex = oldIniIdx;
////                                        recentIdx = index;
//                                    }else{
//                                        iniIndex = index;
//                                        oldIniIdx = index;
////                                        recentIdx = index;
//                                    }
//                                }
//                               
//                                for(int j=0;j<pos2.length;j++){
//                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
//                                    long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                    //** Line 1202 got strange result it shoul produce 24477283769
//                                    if(alnMerMap.containsKey(indexAlnCode)){
//                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode, merList);
//                                        
//                                    }else{
//                                        ArrayList<Long> merList = new ArrayList();
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode,merList);                                        
//                                    }   
//                                }
//                                
//                                /****************************/
//                                
//                                recentIdx = index;
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
                        ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }else{
                        ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
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
                Enumeration<ShortgunSequence> seqsComp = input.getInputSequence().elements();
                while(seqsComp.hasMoreElements()){
                    ShortgunSequence seq = seqsComp.nextElement();                                  // get ShortgunSequence from InputSequence
                    this.alnMerMap = new LinkedHashMap();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
                    String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  
//                    System.out.println("******Input Sequence check " + compSeq);
//                    System.out.print(chr.getName()+" - strand\t"); 
                    
                    Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex> 
                    /******** New Part (fixed wrong mer count) **********/
                    long oldIniIdx = 0;
                    long newIniIdx = 0;
                    long iniIndex = 0;
                    long recentIdx = 0;
                    boolean firstMatchCheck = false;
                    /****************/
                    for(int i=0;i<(compSeq.length()-mer)+1;i++){                                    // Windowing   
                        int index = i;                                                              // index at aligncompliment and non compliment is not different. It not effect any thing. we just know strand notation is enough
                        String sub = compSeq.substring(i, i+mer);
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){
                            m = m<<28;
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2ComplimentV2(m);                            // Do alignment by alignment function specific for compliment sequence
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                            /*------------------------------- add iniIndex in front of 29 bit [strand|position] -----------------------------*/
                            if(pos2 != null){
//                                if(pos2.length == 1){           // (Not work) already check for repeat in same chromosome by checking alignment result must have one match result in this chromosome
//                                for(int j=0;j<pos2.length;j++){
//                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
//                                
//                                       if(alnCodeCheckList.containsKey(alnCode)){
//                                        long iniIndex = alnCodeCheckList.get(alnCode);
//                                        
//                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                        
//                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode, merList);
//                                        
//                                    }else{
//                                        long iniIndex = index;
//                                        alnCodeCheckList.put(alnCode, iniIndex);
//                                        
//                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                        
//                                        ArrayList<Long> merList = new ArrayList();
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode,merList);
//                                        
//                                    }
//                                }
                                

                                /******** New Part (fixed wrong mer count) Version 3 **********/
                                for(int j=0;j<pos2.length;j++){
                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
                                    
                                    
                                    if(alnCodeCheckList.containsKey(alnCode)){
 
                                        iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
                                        
                                        /**
                                         * Case check to solve the problem. In case, when position-index is the same value but actually it different peak.
                                         * To check continuity of this alnCode. We reserve index 0 of merList to store the recent index.
                                         * Check continuity of index from different between recent index and current index.
                                         */
                                        
                                        if(index-merList.get(0)==1){                                // Case check to solve the problem. In case, when position-index is the same value but actually it different peak
                                            iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                            indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                            merList.remove(0);
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode, merList);
                                        }else{
                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode,merList);
                                        }
                                        /*******************/
                                        
                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(0,(long)index);
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode,merList);
                                        
                                    }
    //                                    merList = null;
    //                                    System.gc();                                       
                                }
                                
                                /***************************************************************/


                                /******** New Part (fixed wrong mer count) **********/
                                
                                /* Specify iniIdx check from continueously of index */
//                                newIniIdx = index;
//                                if(firstMatchCheck==false){
//                                    // it's first time
//                                    iniIndex = index;
//                                    oldIniIdx = index;
//                                    firstMatchCheck = true;
//                                }else{
//                                    if(newIniIdx-recentIdx==1){
//                                        iniIndex = oldIniIdx;
//                                    }else{
//                                        iniIndex = index;
//                                        oldIniIdx = index;
//                                    }
//                                }
//                               
//                                for(int j=0;j<pos2.length;j++){
//                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
//                                    long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                    //** Line 1202 got strange result it shoul produce 24477283769
//                                    if(alnMerMap.containsKey(indexAlnCode)){
//                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode, merList);
//                                        
//                                    }else{
//                                        ArrayList<Long> merList = new ArrayList();
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode,merList);                                        
//                                    }   
//                                }
//                                
//                                /****************************/
//                                recentIdx = index;
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
                        
                        ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }else{
                        ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
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
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();
                
                
                
            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        alinResult.addMapResult(this.alnRes);
        /* After Alignment we will loop over ShortgunSequence and do the Alignmentcount */
        //AlignmentResultRead alinResult = new AlignmentResultRead();  
//        Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
//        while(seqs.hasMoreElements()){
//                    ShortgunSequence seq = seqs.nextElement();
//                    seq.countAlignmentData();                                       // count alignment data (this function must be call after information is ready)
//                    alinResult.addResult(seq);                                      // add ShortgunSequence into AlignmentResultRead  
//        }
        
        //return res;
        return alinResult;
    }

    public AlignmentResultRead alignV4(ReferenceSequence ref, InputSequence input, int numMer, int threshold) {
        
        this.setReferenceSequence(ref);
        return alignV4(input,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignV4(InputSequence input, int numMer, int threshold){
        /**
         * The different between alignV3 and V4 is we add another map that contain variant information (SNP,Indel)
         * Suitable data structure V3. 
         */
        
        alnRes = new LinkedHashMap();                                     // initialize this hashmap one time per alinment job
        this.varBox = new LinkedHashMap();
        int numUnMapMer = 0;
        int numUnMapBase = 0;
        //AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;
        
        
        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            try {
                
                
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = SequenceUtil.createAllReference(chr, this.mer);   // Create or import all reference [chromosome reference, repeat index, repeat Marker]
                while(seqs.hasMoreElements()){                                              // Loop over ShortgunSequence contain in InputSequence 
                    ShortgunSequence seq = seqs.nextElement();
                    this.alnMerMap = new LinkedHashMap();                                         // initialize this hashmap every time when start new loop of Shortgun Read
                    this.varType = new LinkedHashMap();
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence
                    
//                    System.out.print(chr.getName()+" + strand\t");           
                    
                    Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
                    /******** New Part (fixed wrong mer count) **********/
                    long oldIniIdx = 0;
                    long newIniIdx = 0;
                    long iniIndex = 0;
                    long recentIdx = 0;
                    boolean firstMatchCheck = false;
                    /****************/
                    for(int i=0;i<(s.length()-mer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]    
                        int index = i;
                        String sub = s.substring(i, i+mer);                                 // cut String sequence into sub string sequence (mer length long) 
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){                                                          
                            m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose 
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2(m);                                // Do alignment with binary search (pos2[] cantain 64 bit long [mer code | position])
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                            /*------------------------------- add iniIndex in front of 29 bit [strand|position] -----------------------------*/ 
                            if(pos2 != null){
                                //if(pos2.length == 1){                // (Not work) already check for repeat in same chromosome by checking alignment result must have one match result in this chromosome
                                
                                
                                /* Old Part *//*
                                for(int j=0;j<pos2.length;j++){
                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
                                    
                                    
                                    if(alnCodeCheckList.containsKey(alnCode)){
                                        long iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode, merList);
                                        
                                    }else{
                                        long iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode,merList);
                                        
                                    }
    //                                    merList = null;
    //                                    System.gc();                                       
                                }
                                */
                                
                                /******** New Part (fixed wrong mer count) Version 3 **********/
                                for(int j=0;j<pos2.length;j++){
                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
                                    
                                    
                                    if(alnCodeCheckList.containsKey(alnCode)){
 
                                        iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
                                        
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
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode, merList);
                                        }else{
                                            /**
                                             * it's not continue. So, iniIndex has change to present index                                                                                  
                                             */
                                            
                                            
//                                            /**
//                                             * Classify for variation type
//                                             * SNP = 1
//                                             * idel = 2
//                                             * other = 3
//                                             */
//                                            numUnMapMer = index - merList.size();
//                                            numUnMapBase = (numUnMapMer-this.mer)+1;
//                                            
//                                            if(numUnMapBase == 0){
//                                                
//                                            }else if{
//                                                
//                                            }hhr
                                            
                                            
                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode,merList);
                                        }
                                        /*******************/
                                        
                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(0,(long)index);
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode,merList);
                                        
//                                        /**
//                                         * Set default variation type
//                                         *  normal = 0
//                                         */
//                                        this.varType.put(indexAlnCode, (byte)0);
 
                                    }
    //                                    merList = null;
    //                                    System.gc();                                       
                                }
                                
                                /***************************************************************/
                                
                                
                                
                                /******** New Part (fixed wrong mer count) Version 2 (not work)**********/
                                
                                /* Specify iniIdx check from continueously of index */
//                                newIniIdx = index;
//                                if(firstMatchCheck==false){
//                                    // it's first time
//                                    iniIndex = index;
//                                    oldIniIdx = index;
////                                    recentIdx = index;
//                                    firstMatchCheck = true;
//                                }else{
//                                    if(newIniIdx-recentIdx==1){
//                                        iniIndex = oldIniIdx;
////                                        recentIdx = index;
//                                    }else{
//                                        iniIndex = index;
//                                        oldIniIdx = index;
////                                        recentIdx = index;
//                                    }
//                                }
//                               
//                                for(int j=0;j<pos2.length;j++){
//                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
//                                    long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                    //** Line 1202 got strange result it shoul produce 24477283769
//                                    if(alnMerMap.containsKey(indexAlnCode)){
//                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode, merList);
//                                        
//                                    }else{
//                                        ArrayList<Long> merList = new ArrayList();
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode,merList);                                        
//                                    }   
//                                }
//                                
//                                /****************************/
//                                
//                                recentIdx = index;
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
                        ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }else{
                        ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
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
                Enumeration<ShortgunSequence> seqsComp = input.getInputSequence().elements();
                while(seqsComp.hasMoreElements()){
                    ShortgunSequence seq = seqsComp.nextElement();                                  // get ShortgunSequence from InputSequence
                    this.alnMerMap = new LinkedHashMap();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
                    String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  
//                    System.out.println("******Input Sequence check " + compSeq);
//                    System.out.print(chr.getName()+" - strand\t"); 
                    
                    Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex> 
                    /******** New Part (fixed wrong mer count) **********/
                    long oldIniIdx = 0;
                    long newIniIdx = 0;
                    long iniIndex = 0;
                    long recentIdx = 0;
                    boolean firstMatchCheck = false;
                    /****************/
                    for(int i=0;i<(compSeq.length()-mer)+1;i++){                                    // Windowing   
                        int index = i;                                                              // index at aligncompliment and non compliment is not different. It not effect any thing. we just know strand notation is enough
                        String sub = compSeq.substring(i, i+mer);
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){
                            m = m<<28;
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2ComplimentV2(m);                            // Do alignment by alignment function specific for compliment sequence
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
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                            /*------------------------------- add iniIndex in front of 29 bit [strand|position] -----------------------------*/
                            if(pos2 != null){
//                                if(pos2.length == 1){           // (Not work) already check for repeat in same chromosome by checking alignment result must have one match result in this chromosome
//                                for(int j=0;j<pos2.length;j++){
//                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
//                                
//                                       if(alnCodeCheckList.containsKey(alnCode)){
//                                        long iniIndex = alnCodeCheckList.get(alnCode);
//                                        
//                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                        
//                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode, merList);
//                                        
//                                    }else{
//                                        long iniIndex = index;
//                                        alnCodeCheckList.put(alnCode, iniIndex);
//                                        
//                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                        
//                                        ArrayList<Long> merList = new ArrayList();
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode,merList);
//                                        
//                                    }
//                                }
                                

                                /******** New Part (fixed wrong mer count) Version 3 **********/
                                for(int j=0;j<pos2.length;j++){
                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
                                    
                                    
                                    if(alnCodeCheckList.containsKey(alnCode)){
 
                                        iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
                                        
                                        /**
                                         * Case check to solve the problem. In case, when position-index is the same value but actually it different peak.
                                         * To check continuity of this alnCode. We reserve index 0 of merList to store the recent index.
                                         * Check continuity of index from different between recent index and current index.
                                         */
                                        
                                        if(index-merList.get(0)==1){                                // Case check to solve the problem. In case, when position-index is the same value but actually it different peak
                                            iniIndex = alnCodeCheckList.get(alnCode);
                                        
                                            indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                            merList.remove(0);
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode, merList);
                                        }else{
                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode,merList);
                                        }
                                        /*******************/
                                        
                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);
                                        
                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                        
                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(0,(long)index);
                                        merList.add(m);
                                        this.alnMerMap.put(indexAlnCode,merList);
                                        
                                    }
    //                                    merList = null;
    //                                    System.gc();                                       
                                }
                                
                                /***************************************************************/


                                /******** New Part (fixed wrong mer count) **********/
                                
                                /* Specify iniIdx check from continueously of index */
//                                newIniIdx = index;
//                                if(firstMatchCheck==false){
//                                    // it's first time
//                                    iniIndex = index;
//                                    oldIniIdx = index;
//                                    firstMatchCheck = true;
//                                }else{
//                                    if(newIniIdx-recentIdx==1){
//                                        iniIndex = oldIniIdx;
//                                    }else{
//                                        iniIndex = index;
//                                        oldIniIdx = index;
//                                    }
//                                }
//                               
//                                for(int j=0;j<pos2.length;j++){
//                                    long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]
//                                    long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
//                                    //** Line 1202 got strange result it shoul produce 24477283769
//                                    if(alnMerMap.containsKey(indexAlnCode)){
//                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode, merList);
//                                        
//                                    }else{
//                                        ArrayList<Long> merList = new ArrayList();
//                                        merList.add(m);
//                                        this.alnMerMap.put(indexAlnCode,merList);                                        
//                                    }   
//                                }
//                                
//                                /****************************/
//                                recentIdx = index;
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
                        
                        ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }else{
                        ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
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
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();
                
                
                
            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        alinResult.addMapResult(this.alnRes);
        /* After Alignment we will loop over ShortgunSequence and do the Alignmentcount */
        //AlignmentResultRead alinResult = new AlignmentResultRead();  
//        Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
//        while(seqs.hasMoreElements()){
//                    ShortgunSequence seq = seqs.nextElement();
//                    seq.countAlignmentData();                                       // count alignment data (this function must be call after information is ready)
//                    alinResult.addResult(seq);                                      // add ShortgunSequence into AlignmentResultRead  
//        }
        
        //return res;
        return alinResult;
    }
    
    public AlignmentResultRead alignV5(ReferenceSequence ref, InputSequence input, int numMer, int threshold) {
        
        this.setReferenceSequence(ref);
        return alignV3(input,numMer,threshold);
        
    }
    
    public AlignmentResultRead alignV5(InputSequence input, int numMer, int threshold){
        /**
         * The different between alignV5 and V3 is we use repeat marker to align repeat
         */
        alnRes = new LinkedHashMap();                                     // initialize this hashmap one time per alinment job
       
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        this.mer = numMer;
        
        
        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){                                                      // Loop chromosome contain in ReferenceSequence
            try {
                
                
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = SequenceUtil.createAllReferenceV2(chr, this.mer, 'a');   // Create or import all reference [chromosome reference, repeat index, repeat Marker]               
                while(seqs.hasMoreElements()){                                              // Loop over ShortgunSequence contain in InputSequence 
                    Map<Integer,ArrayList<Integer>> linkIndexCheck = new LinkedHashMap();                       // HashMap contain data that has been use to check for repeat jump
                    ShortgunSequence seq = seqs.nextElement();
                    this.alnMerMap = new LinkedHashMap();                                         // initialize this hashmap every time when start new loop of Shortgun Read

                    String s = seq.getSequence();                                           // get String sequence of selected ShortgunSequence

                    Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex>
                    /******** New Part (fixed wrong mer count) **********/
                    long oldIniIdx = 0;
                    long newIniIdx = 0;
                    long iniIndex = 0;
                    long recentIdx = 0;
                    boolean firstMatchCheck = false;
                    /****************/
                    for(int i=0;i<(s.length()-mer)+1;i++){                                  // Big window (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]    
                        int index = i;
                        String sub = s.substring(i, i+mer);                                 // cut String sequence into sub string sequence (mer length long) 
                        long m = SequenceUtil.encodeMer(sub, mer);                          // encode sub string sequence (code is 36 bit max preserve the rest 28 bit for position)
                        if(m!=-1){                                                          
                            m = m<<28;                                                      // shift left 28 bit for optimization binary search purpose
                            
                            /**
                             * write new implement code here
                             * align repeat first
                             * if output is nothing goes to normal align
                             * align repeat should pass s (a sequence string) and all currently move information because we have to move small window inside it
                             */
                            
                            long posR[] = encoded.align3(m, s, index, numMer, linkIndexCheck);
//                            ArrayList<Long> posR = encoded.align3(m, s, index, numMer, linkIndexCheck);
                            if(posR == null){

                                long pos2[] = encoded.align2(m);                                // Do alignment with binary search (pos2[] cantain 64 bit long [mer code | position])
                                long pos = -1;
                                if(pos2!=null&&pos2.length>0){
                                    pos = pos2[0];
                                    pos = pos2.length;
                                }

                                int idx = (int) (pos-i);
                                if(pos<0){
                                  idx = 0;
                                }

                                int totalMer = (seq.getShortgunLength()-mer)+1;

                                /*************************************************************************************************************/
                                /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                                /*------------------------------- add iniIndex in front of 29 bit [strand|position] -----------------------------*/ 
                                if(pos2 != null){

                                    /******** New Part (fixed wrong mer count) Version 3 **********/
                                    for(int j=0;j<pos2.length;j++){
                                        long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]


                                        if(alnCodeCheckList.containsKey(alnCode)){

                                            iniIndex = alnCodeCheckList.get(alnCode);

                                            long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);

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
                                                merList.add(0,(long)index);
                                                merList.add(m);
                                                this.alnMerMap.put(indexAlnCode, merList);
                                            }else{
                                                /**
                                                 * it's not continue. So, iniIndex has change to present index                                                                                  
                                                 */

                                                iniIndex = index;
                                                alnCodeCheckList.put(alnCode, iniIndex);

                                                indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                                merList = new ArrayList();
                                                merList.add(0,(long)index);
                                                merList.add(m);
                                                this.alnMerMap.put(indexAlnCode,merList);
                                            }
                                            /**************************************************************************************************/

                                        }else{
                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            ArrayList<Long> merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode,merList);

                                        }

                                    }
                                    /***************************************************************/
                                }                         
                                /*************************************************************************************************************/
                            }else if(posR!=null){
                                /******** New Part (fixed wrong mer count) Version 3 **********/
                                long mask29Bit = 536870911;
                                for(int j=0;j<posR.length;j++){
                                    if(posR[j] == 0){
                                        /**
                                         * Check for 0 value. posR Array store ~39bit of [merCount|strand|position] if it has value at least merCount must = 1. So, it impossible to have zero element.
                                         * Also, the posR has design to keep only the element that have information and store at first element first then so on by order. 
                                         * This mean if we found the first element that have 0 value the element follow by this is all 0 as well. It useless to continue looping. 
                                         * So, we can break the loop to reduce computational time. 
                                         */

                                        break;
                                    }
                                    
                                    int merCount = (int)(posR[j]>>29);       
                                    long alnCode = (posR[j]&mask29Bit) - index;     // posR is ~39 bit [merCount|strand|position] ; algncode is 29 bit [strand|alignPosition]. alignposition is position - index


                                    if(alnCodeCheckList.containsKey(alnCode)){

                                        iniIndex = alnCodeCheckList.get(alnCode);

                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);

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
                                            merList.add(0,(long)index);
                                            for(int num=0;num<merCount;num++){
                                                merList.add(0L);
                                            }
                                            
                                            this.alnMerMap.put(indexAlnCode, merList);
                                        }else{
                                            /**
                                             * it's not continue. So, iniIndex has change to present index                                                                                  
                                             */

                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            for(int num=0;num<merCount;num++){
                                                merList.add(0L);
                                            }
                                            this.alnMerMap.put(indexAlnCode,merList);
                                        }
                                        /**************************************************************************************************/

                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);

                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(0,(long)index);
                                        for(int num=0;num<merCount;num++){
                                            merList.add(0L);
                                        }
                                        this.alnMerMap.put(indexAlnCode,merList);

                                    }

                                }
                            }
                            
                        } 
                                            
                    }
                    
                    /*************************************************************************************************************/
                    /* -------------------------New Implement Part Cons. (Not Stroe in object)---------------------------------------------*/
                   
                    if(this.alnRes.containsKey(seq.getReadName())){                     // Check for existing of ReadName (if exist put result code on existing ArrayList<Long>
                        ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }else{
                        ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 29 bit compose of [strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add count number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }
                    /*-----------------------------------------------------------------------------------------------------------*/
                    /*************************************************************************************************************/
                    

                }
                
                /*-------------------- Do compliment alignment -------------------------------*/
                /* Do the same algorithm but use function for compliment */
                Enumeration<ShortgunSequence> seqsComp = input.getInputSequence().elements();
                while(seqsComp.hasMoreElements()){
                    Map<Integer,ArrayList<Integer>> linkIndexCheck = new LinkedHashMap();                       // HashMap contain data that has been use to check for repeat jump
                    ShortgunSequence seq = seqsComp.nextElement();                                  // get ShortgunSequence from InputSequence
                    this.alnMerMap = new LinkedHashMap();

                    String s = seq.getSequence();                                                   // get sequence form ShortgunSequence
                    String invSeq = SequenceUtil.inverseSequence(s);                                // Do invert sequence (ATCG => GCTA)
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);                       // Do compliment on invert sequence (GCTA => CGAT)  

                    Map<Long,Long> alnCodeCheckList = new HashMap();                    // This map is a checklist for alncode to indicate the iniIndex Map<Long,Long> => Map<alnCode,iniIndex> 
                    
                    
                    
                    
                    
                    /******** New Part (fixed wrong mer count) **********/
                    long oldIniIdx = 0;
                    long newIniIdx = 0;
                    long iniIndex = 0;
                    long recentIdx = 0;
                    boolean firstMatchCheck = false;
                    /****************/
                    for(int i=0;i<(compSeq.length()-mer)+1;i++){                                    // Windowing   
                        int index = i;                                                              // index at aligncompliment and non compliment is not different. It not effect any thing. we just know strand notation is enough
                        String sub = compSeq.substring(i, i+mer);
                        long m = SequenceUtil.encodeMer(sub, mer);
                        if(m!=-1){
                            m = m<<28;
                            long posR[] = encoded.align3Compliment(m, compSeq, index, numMer, linkIndexCheck);
//                            ArrayList<Long> posR = encoded.align3Compliment(m, compSeq, index, numMer, linkIndexCheck);

                            if(posR == null){
                            
                                long pos2[] = encoded.align2ComplimentV2(m);                            // Do alignment by alignment function specific for compliment sequence
                                long pos = -1;
                                if(pos2!=null&&pos2.length>0){
                                    pos = pos2[0];
                                    pos = pos2.length;
                                }
                            
                                int idx = (int) (pos-i);
                                if(pos<0){
                                  idx = 0;
                                }
                                int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                                /*************************************************************************************************************/
                                /* -------------------------New Implement Part (Not Stroe in object)---------------------------------------------*/
                                /*------------------------------- add iniIndex in front of 29 bit [strand|position] -----------------------------*/
                                if(pos2 != null){
                                    /******** New Part (fixed wrong mer count) Version 3 **********/
                                    for(int j=0;j<pos2.length;j++){
                                        long alnCode = pos2[j] - index;     // pos is 29 bit [strand|position] ; algncode is 29 bit [strand|alignPosition]


                                        if(alnCodeCheckList.containsKey(alnCode)){

                                            iniIndex = alnCodeCheckList.get(alnCode);

                                            long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);

                                            /**
                                             * Case check to solve the problem. In case, when position-index is the same value but actually it different peak.
                                             * To check continuity of this alnCode. We reserve index 0 of merList to store the recent index.
                                             * Check continuity of index from different between recent index and current index.
                                             */

                                            if(index-merList.get(0)==1){                                // Case check to solve the problem. In case, when position-index is the same value but actually it different peak
                                                iniIndex = alnCodeCheckList.get(alnCode);

                                                indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)
                                                merList.remove(0);
                                                merList.add(0,(long)index);
                                                merList.add(m);
                                                this.alnMerMap.put(indexAlnCode, merList);
                                            }else{
                                                iniIndex = index;
                                                alnCodeCheckList.put(alnCode, iniIndex);

                                                indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                                merList = new ArrayList();
                                                merList.add(0,(long)index);
                                                merList.add(m);
                                                this.alnMerMap.put(indexAlnCode,merList);
                                            }
                                            /*******************/

                                        }else{
                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            ArrayList<Long> merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            merList.add(m);
                                            this.alnMerMap.put(indexAlnCode,merList);

                                        }

                                    }

                                    /***************************************************************/

                                }
                            }else if(posR != null){
                                /******** New Part (fixed wrong mer count) Version 3 **********/
                                long mask29Bit = 536870911;
                                for(int j=0;j<posR.length;j++){
                                    if(posR[j] == 0){
                                        /**
                                         * Check for 0 value. posR Array store ~39bit of [merCount|strand|position] if it has value at least merCount must = 1. So, it impossible to have zero element.
                                         * Also, the posR has design to keep only the element that have information and store at first element first then so on by order. 
                                         * This mean if we found the first element that have 0 value the element follow by this is all 0 as well. It useless to continue looping. 
                                         * So, we can break the loop to reduce computational time. 
                                         */

                                        break;
                                    }
                                    
                                    
                                    
                                    int merCount = (int)(posR[j]>>29);       
                                    long alnCode = (posR[j]&mask29Bit) - index;     // posR is ~39 bit [merCount|strand|position] ; algncode is 29 bit [strand|alignPosition]. alignposition is position - index


                                    if(alnCodeCheckList.containsKey(alnCode)){

                                        iniIndex = alnCodeCheckList.get(alnCode);

                                        long indexAlnCode = (iniIndex<<29)+alnCode;                 // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                        ArrayList<Long> merList = this.alnMerMap.get(indexAlnCode);

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
                                            merList.add(0,(long)index);
                                            for(int num=0;num<merCount;num++){
                                                merList.add(0L);
                                            }
                                            
                                            this.alnMerMap.put(indexAlnCode, merList);
                                        }else{
                                            /**
                                             * it's not continue. So, iniIndex has change to present index                                                                                  
                                             */

                                            iniIndex = index;
                                            alnCodeCheckList.put(alnCode, iniIndex);

                                            indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                            merList = new ArrayList();
                                            merList.add(0,(long)index);
                                            for(int num=0;num<merCount;num++){
                                                merList.add(0L);
                                            }
                                            this.alnMerMap.put(indexAlnCode,merList);
                                        }
                                        /**************************************************************************************************/

                                    }else{
                                        iniIndex = index;
                                        alnCodeCheckList.put(alnCode, iniIndex);

                                        long indexAlnCode = (iniIndex<<29)+alnCode;                  // indexAlnCode has 37 bit [iniIndex|Strand|Position] iniIndex(8bit),Strnd(1bit),Position(28bit)

                                        ArrayList<Long> merList = new ArrayList();
                                        merList.add(0,(long)index);
                                        for(int num=0;num<merCount;num++){
                                            merList.add(0L);
                                        }
                                        this.alnMerMap.put(indexAlnCode,merList);

                                    }
                                }
                            }
                                                   
                            /*-----------------------------------------------------------------------------------------------------------*/
                            /*************************************************************************************************************/
                        }
                                           
                    }
                    
                    /*************************************************************************************************************/
                    /* -------------------------New Implement Part Cons. (Not Stroe in object)---------------------------------------------*/
                   
                    if(this.alnRes.containsKey(seq.getReadName())){                     // Check for existing of ReadName (if exist put result code on existing ArrayList<Long>
                        
                        ArrayList<Long> countChrIdxStrandAlnList = this.alnRes.get(seq.getReadName()); //get existing Arraylist
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }else{
                        ArrayList<Long> countChrIdxStrandAlnList = new ArrayList();
                        Set keySet = this.alnMerMap.keySet();
                        Iterator keyIter =keySet.iterator();
                        while(keyIter.hasNext()){
                            long idxStrandAln = (long)keyIter.next();                      // strandAln has 37 bit compose of [iniIndex|strand|alignPosition]
                            long count = this.alnMerMap.get(idxStrandAln).size()-1;          // we can get number of count from number of member in merList and should minus with 1 (because index 0 has been reseve for checking index continuity)
                            long chrIdxStrandAln = (chr.getChrNumber()<<37)+idxStrandAln;     // shift left 37 bit beacause we want to add chr number on the front of strandAln which has 37 bit
                            long countChrIdxStrandAln = (count<<42)+chrIdxStrandAln;          // shift left 42 bit beacause we want to add count number on the front of chrStrandAln which has 42 bit 
                            
                            if(count>threshold){                                                // case check to filter small count peak out (use user specify threshold)
                                countChrIdxStrandAlnList.add(countChrIdxStrandAln);
                            }
                        }
                        this.alnRes.put(seq.getReadName(), countChrIdxStrandAlnList);
                    }

                    /*-----------------------------------------------------------------------------------------------------------*/
                    /*************************************************************************************************************/

                }
                /* End */
                
                encoded.lazyLoad();         // clear memmory
                encoded = null;
                
                System.gc();
                
                
                
            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        alinResult.addMapResult(this.alnRes);
        
        return alinResult;
    }
}