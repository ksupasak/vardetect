/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import biotec.bsi.ngs.vardetect.core.CombineReferenceSequence.FixedDuplicateMer;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeMer;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import static htsjdk.samtools.util.SequenceUtil.reverseComplement;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author soup
 */
public class CombineReferenceSequence extends ReferenceSequence  implements ThreadPoolCallBack {
    
    
    boolean debug = false;
    String output = null;
    String foutput = null;
    PrintStream pout =null;
    PrintStream fout =null;
    
    long skip_read = 0;
    long total_read = 0;
    
    int number_of_thread = 1;
    
    int chunk_read = 0;
    
    int minimum_peak_pattern1 = 10;
    int minimum_peak_pattern2 = 5;
    
    int filter_mer_coverage = 50;
    int filter_mode = 1; // junction_only = 1, indel = 2, junction+indel = 3
    int filter_ref_repeat = 50;
    int filter_dup_count = 1;
    
    int minimumIndelSize = 1;     // minimum size of indel Eg 1 mean allow indel has indel size more or equal to 2

    public void setMinimumIndelSize(int minimumIndelSize) {
        this.minimumIndelSize = minimumIndelSize;
    }
    
    public void setFilterMerCoverage(int f){
        this.filter_mer_coverage = f;
    }
    
    public void setFilterMode(int f){
        this.filter_mode = f;
    }
    
    public void setFilterRefRepeatCount(int f){
        this.filter_ref_repeat = f;
    }
    
    public void setFilterDupCount(int f){
        this.filter_dup_count = f;
    }
    
    public void setChunkRead(int c){
        this.chunk_read = c;
    }
    
    public void setNumberOfThread(int t){
        number_of_thread = t;
    }
    public void setMinimumPeakPattern(int min1, int min2){
        minimum_peak_pattern1 = min1;
        minimum_peak_pattern2 = min2;
    }
    public void setSkipRead(long s){
        skip_read = s;
    }
    public void setTotalRead(long s){
        total_read = s;
    }
    
    public void setOutputFile(String s){
        this.output = s;
    }
    
    public void setOutputSVFile(String s){
        this.foutput = s;
    }
    
    
    ThreadPool pool = null;
    HashMap queue = new HashMap();
    
    
    public void finishBatch(int i){
        queue.remove(i);
        if(queue.size()==0)pool.shutdown();
    }
    
    public void runProfileSV(String seq_file) throws FileNotFoundException{
        
        File bamFile = new File(seq_file); 
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader samReader = samReaderFactory.open(bamFile);

        pout = System.out;
                
        if(output!=null){
            
            pout = new PrintStream(new FileOutputStream(output));
            
        }
        fout = System.out;
        
        if(foutput!=null){
            fout = new PrintStream(new FileOutputStream(foutput));
        }
        
        boolean finished=false;
        
        if(this.total_read>0){
            
        pool = new ThreadPool(number_of_thread,this);
        long time = this.total_read/number_of_thread;
        
        for (int i = 0; i < number_of_thread ; i++) {
            SVProfiler task = new SVProfiler(this,seq_file ,time,time*i+skip_read, i);
            queue.put(i, task);
            pool.execute(task, i);
        }}
        else{
            pool = new ThreadPool(number_of_thread,this);

           SAMRecordIterator it = samReader.iterator();
           
           int count = 0 ;
           int i = 0 ;
           long time = this.chunk_read;
           while(it.hasNext()){
              it.next();
              if(count!=0&&count%this.chunk_read==0){
                  System.out.println("Chunk "+i+" "+time*i);
                  SVProfiler task = new SVProfiler(this,seq_file ,time,time*i+skip_read, i);
                  queue.put(i, task);
                  pool.execute(task, i);

                  i++;
              }
              count++; 
           }
            SVProfiler task = new SVProfiler(this,seq_file ,count%this.chunk_read-1,time*i+skip_read, i);
            queue.put(i, task);
            pool.execute(task, i);            
               
        }
        
        while(true){
            if(pool.isAllJobDone()){
                System.out.println("All job Done");
                break;
            }
        }
   
    }

    public void prepare() throws IOException {
        if(this.random_access==false){
            
            loadAllSequence();
            
            
        }
        
        searchMer(0);

    }

    public void loadAllSequence() {

       /**
             * Extract chromosome and create index file
             */
            this.chrs.clear();
            
            Charset charset = Charset.forName("US-ASCII");
            Path path = Paths.get(filename);
            String chr = null;

            StringBuffer seq = new StringBuffer();

            try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
                String line = null;

                while ((line = reader.readLine()) != null) {

                    if(line.charAt(0)=='>'){

                        if(chr!=null){

                            System.out.println("CfHR : "+chr+" Size : "+seq.length());
                              
                            ChromosomeSequence c = new ChromosomeSequence(this,chr,seq);
                            this.addChromosomeSequence(c);
                            
                            
                        }
                        seq = new StringBuffer();
                        chr = line.substring(1,line.length());
                       chr = chr.split(" ")[0];
                     
                    }
                     else{
                        seq.append(line.trim());
                    }
                }
                
                 if(seq.length()>0){

                    System.out.println("CHR : "+chr+" Size : "+seq.length());
                    ChromosomeSequence c = new ChromosomeSequence(this,chr,seq);

                    seq = null;

                    this.addChromosomeSequence(c);
                }

                
            }catch (IOException x) {
                System.err.format("IOException: %s%n", x);
            } 


    }
    
    
    public class SVProfiler implements Runnable {

    long skip_read = 0;
    long total_read = 0;   
    String seq_file = null;
    int thread_id = 0 ;
    CombineReferenceSequence ref;
    
    public SVProfiler(CombineReferenceSequence ref,String seq_file, long count,long skip,int thread_id){
        skip_read = skip;
        total_read = count;
        this.seq_file = seq_file;
        this.ref = ref;
        this.thread_id = thread_id;
        
    }  
        
    public void run() {
        try {
            ref.profileSV(seq_file, total_read, skip_read, thread_id);
        } catch (IOException ex) {
            Logger.getLogger(CombineReferenceSequence.class.getName()).log(Level.SEVERE, null, ex);
        }
    }


    }
    public static String reverseComplement(String s){
        return htsjdk.samtools.util.SequenceUtil.reverseComplement(s);
    }
    
    class MappedResult{
        int direction = 0;
        int[] pos;
        
        MappedResult(int[] pos, int direction){
            this.pos = pos;
            this.direction = direction;
        }
        
        public int[] getPos(){
            return pos;
        }
       
        public int getDirection(){
            return this.direction;
        }
    }
    
    public  void profileSV(String filename,long countr,long skip, int thread_id) throws IOException{
        
        
        
        
        File bamFile = new File(filename); 
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader samReader = samReaderFactory.open(bamFile);

        
        
        SAMRecordIterator rr = samReader.iterator();
        int last_percent = 0 ;
        for(long z=0;z<countr+skip-1;z++){
            
            if(z-skip>0){
            int percent = (int)((z-skip)*100/countr);
            
            if(percent>0&&percent!=last_percent){        
                if(this.random_access==false&&this.output!=null)System.out.println("Process "+thread_id+ " "+percent+"%");
                last_percent = percent;
            }}
            
            SAMRecord r = rr.next();
            if(z<skip-1)continue;
            try{
                
            String read = r.getReadString();
            /**
             * check read unmap or not 
             * if read is unmap => flagUnmap will be 1
             * else read is map => flagUnmap will be 0
             */
            byte flagUnmap = 1;
            if(r.getReadUnmappedFlag()==true){
                flagUnmap = 1;
            }else{
                flagUnmap = 0;
            }
            /*******************************************/        
                    
            String sb = read;
            
              
            int kmer = 16;
            int sliding = 1;
            
            int n = (sb.length()-kmer)/sliding;       
            long cmer = -1;
            long mask = 0; 
            int count =0;
            int start =0;

            HashMap rx = new HashMap();
            
            
            String max="";
            String max2="";
            
            
            
            ArrayList res = new ArrayList(r.getReadLength()*2);
            
            
            
            for(int i =0;i<kmer;i++)mask=mask*4+3;

            for(int di = 0;di<2;di++){
            
            if(di==1){
                sb = reverseComplement(sb);
                
            }    
                
            
            for(int i =0;i<n;i++){

                long pos = i;
                
                
                char chx = sb.charAt(i*sliding+kmer-1);
                if(chx!='N'){
                    if(cmer==-1){
                        String s = sb.substring(i*sliding,i*sliding+kmer);
                        cmer = encodeMer(s,kmer);
                    }else{

                        int t =-1;
                        switch(chx){
                            case 'A':
                            case 'a':
                                t=0; // 00
                                break;
                            case 'T':
                            case 't': 
                                t=3; // 11
                                break;
                            case 'U':
                            case 'u': 
                                t=3; // 11 (in case of RNA T has change to U)
                                break;
                            case 'C':
                            case 'c':
                                t=1; // 01 
                                break;
                            case 'G':
                            case 'g':
                                t=2; // 10 
                                break;
                            default : 
                                t=-1;
                            break;

                        }
                        if(t>=0){

                            cmer *= 4;
                            cmer &= mask;
                            cmer += t;

                        }else{
                            cmer = -1;
                            i+=kmer;
                        }  
                    }  
                    
                    
//                    if(i%10000000==0)System.out.println("Encode - "+name+" "+(i*sliding));

                    if(true||cmer>=0){  
//                        long x = (cmer<<(64- kmer*2))|pos;
//                        long x = (cmer<<32)|pos;        // left Shift 28 bit is mean we reserve 28 bit for position (we fixed it) 
//                        list[i] = x;
                        int[] result = this.searchMer((int)cmer,thread_id);
                        
//                        res.add(result);
                        
                        res.add(new MappedResult(result,di));
                        
                        
                        if(result!=null){
//                            System.out.print(""+i+ " "+result.length);
                            
                            for(int j=0;j<result.length;j++){
                            
                                pos = (((long)result[j])&mask);
                                
                                if(debug)System.out.print(" "+di+"-"+pos);
                                
                                
                                long xpos = pos - i;
                                String key = ""+xpos+":"+di;
                                
                                if(rx.containsKey(key)){
                                    int v = (int)rx.get(key)+1;
                                    rx.put(key,v);
                                    if(rx.get(max)==null)max=key;
                                    if(v>(int)rx.get(max)){
                                        max = key;
                                    }
                                           
                                    
                                }else
                                   rx.put(key, 1);
                                
                            }
                            
                             if(debug)System.out.println(" "+di);
                            
//                            for(int j=0;j<result.length;j++)  System.out.print(" "+result[j]);
                            
//                            System.out.println();
                            
//                            System.out.print(result.length);
                          
                            
                        }else{
                             if(debug)System.out.println(" "+di+"-");
                        }
                        
                        

                        count++;
                    }
                    
                

                }
                
              
                 
            }
            
                
            }
            
        //            rx.entrySet()
            
//            Iterator it = rx.entrySet().iterator();
////            m
//            while (it.hasNext()) {
//                    
//                    Map.Entry pair = (Map.Entry)it.next();
//                    String k = (String)pair.getKey();
//                    int v = (int)pair.getValue();
//                    if(k.compareTo(max)!=0){
//                        if(rx.get(max2)==null)max2=k;
//                        else{
//                            if(v>(int)rx.get(max2)){
//                                max2=k;
//                            }
//                        }
//                        
//                    }
////                    it.remove(); // avoids a ConcurrentModificationException
//            }
             
            
     
             
            
            if(rx.get(max)!=null&&(int)rx.get(max)>minimum_peak_pattern1){
             
                
             long pmax = Long.valueOf(max.split(":")[0]);
             int dix = Integer.valueOf(max.split(":")[1]);
             
             
             
   
             CombinedPos cp = this.getChromosomePos(pmax);
             
             String chrA = cp.getChr();
             
             int posA = cp.getPos();
             String refseq,refseq2="";
             
             if(!this.random_access){
                refseq = this.getReferenceSequence(cp.getChr(), cp.getPos(), sb.length());
//                refseq2 = this.getReferenceSequence(cp2.getChr(), cp2.getPos(), 150);
                if(refseq == null){
                    continue;
                }

             }else{
                refseq = this.getReferenceSequenceFromIndex(cp.getChr(), cp.getPos(), sb.length());
//                refseq2 = this.getReferenceSequenceFromIndex(cp2.getChr(), cp.getPos(), 150);
                if(refseq == null){
                    continue;
                }
             }
             String x = refseq;
             if(dix==1)
             refseq = reverseComplement(refseq);
             
//             if(dix==1)
//             read = reverseComplement(read);
             

             
             
             
                EvaluateResult resc = evaluateConfirmation(read,refseq,null);
               
           
                
                HashMap sc = new HashMap();
                String rescseq = resc.getEvaluateTag();
                for(int u=0;u<rescseq.length();u++){
                    
                    if(rescseq.charAt(u)=='_'){
                        int p = u;
                        if(dix==1)p=n-p-1;
                        sc.put(p, true);
                    }
                    
                }
 
                HashMap rx2 = new HashMap();
                max2 = "";
                for(int dic=0;dic<2;dic++){
                    
                for(int u=0;u<n;u++){
                    
                    try{
//                     System.out.println(""+(dic*n+j));
                    MappedResult resx = (MappedResult)res.get(dic*n+u);
                    
                    if(resx.getPos()!=null&&sc.containsKey(u)==true){
                        
//                          System.out.println(""+dic+" "+u+" "+" : "+resx.getPos().length);
                            int result[] = resx.getPos();
                            
                            for(int j=0;j<result.length;j++){
                            
                                long pos = (((long)result[j])&mask);
                                long xpos = pos - u;
                                String key = ""+xpos+":"+dic;
//                                 System.out.println("Key : "+key);
                                if(rx2.containsKey(key)){
                                    
                                    int v = (int)rx2.get(key)+1;
                                    
                                    rx2.put(key,v);
                                    if(rx2.get(max2)==null)max2=key;
                                    if(v>(int)rx2.get(max2)){
                                        max2 = key;
                                    }
                                           
                                    
                                }else
                                   rx2.put(key, 1);
                                
                            }
    
                    }
                    }catch(Exception e){
                        break;
//                        System.err.println("xx"+e.getMessage());
//                        
//                        StackTraceElement[] elements = e.getStackTrace();
//                for(int i=0; i<elements.length; i++) {
//                System.err.println(elements[i]);
//                }
                    }
                    
                }
                
                    
                   
                    
                    
                    
                    
                     
                }
                
                 if(rx2.get(max2)!=null&&(int)rx2.get(max2)>minimum_peak_pattern2){

                
                    
//                     System.out.println("Max2 :"+(max2));
                     
                   
                    long pmax2 = Long.valueOf(max2.split(":")[0]);
                    CombinedPos cp2 = this.getChromosomePos(pmax2);
                    String chrB = cp2.getChr();
                    int posB = cp2.getPos();
                    int dix2 = Integer.valueOf(max2.split(":")[1]);

                    
                    if(!this.random_access){
//                        refseq = this.getReferenceSequence(cp.getChr(), cp.getPos(), 150);
                        refseq2 = this.getReferenceSequence(chrB, posB, sb.length());
                        if(refseq2 == null){
                            continue;
                        }
                    }else{
//                        refseq = this.getReferenceSequenceFromIndex(cp.getChr(), cp.getPos(), 150);
                        refseq2 = this.getReferenceSequenceFromIndex(chrB, posB, sb.length());
                        if(refseq2 == null){
                            continue;
                        }
                    }
                    
                    if(dix2==1)
                    refseq2 = reverseComplement(refseq2);
                    
                    //=============================
                    
                    
//                    EvaluateResult resc2 = evaluateConfirmation(read,refseq2,null);
                    EvaluateResult resf = evaluateConfirmation(read,refseq,refseq2);
                    
                    int maxAc = (int)rx.get(max);
                    int maxBc = (int)rx2.get(max2);
                    String chrposA = ""+(chrA+":"+posA);
                    String chrposB = ""+(chrB+":"+posB);
                    String maxposA = max;
                    String maxposB = max2;
                    
                            
                            
                    if(resf.isReversePeak()){
//                        System.out.println("Reverse");
                        resf = evaluateConfirmation(read,refseq2,refseq);
                        
                        String trefseq = refseq;
                        refseq = refseq2;
                        refseq2 = trefseq;
                        
                        int tmaxc = maxAc;
                        maxAc = maxBc;
                        maxBc = tmaxc;
                        
                        String tchrpos = chrposA;
                        chrposA  = chrposB;
                        chrposB = tchrpos;
                        
                        String tmaxpos = maxposA;
                        maxposA = maxposB;
                        maxposB = tmaxpos;
                        
                        
                        
                    }
                    
                    
                    
                    resf.resolve();
                    if(debug){
                        
                        System.out.println("Read        : "+read);
                        System.out.println("Read Length : "+read.length());
                        System.out.println("Mer Length  : "+n);
                        System.out.println("Max1        : "+max);
                        System.out.println("Max1 count  : "+rx.get(max));
                        System.out.println("Max1 fw/rv  : "+dix);
                        System.out.println("Max1 Resc   : "+resc.getEvaluateTag());
                        System.out.println("Read        : "+read);
                        System.out.println("Ref1        : "+refseq);
                        System.out.println("Max2        : "+max2);
                        System.out.println("Max2 count  : "+rx2.get(max2));
                        System.out.println("Max2 fw/rv  : "+dix2);
                        System.out.println("Read        : "+read);
                        System.out.println("Ref2        : "+refseq2);
                        
                  
                        System.out.println("Final       : "+resf.getEvaluateTag());
                        System.out.println("Repeat Count: "+resf.getRefRepeatCount());
                        System.out.println("Resolve     : "+resf.getResolveTag());
                        int bp[] = resf.findBreakPoint();
                        System.out.println("Breakpoint  : "+bp[0]+"-"+bp[1]+" "+read.substring(bp[0], bp[1]) );
                    
                    }
                    String sread = resf.getEvaluateTag();
                       
                    pout.println(""+(z+1)+"\t"+resf.getCoverage()+"\t"+resf.getSwitchCount()+"\t"+resf.getDupCount()+"\t"+maxAc+"\t"+maxBc+"\t"+maxposA+"\t"+maxposB+"\t"+chrposA+"\t"+chrposB+"\t"+r.getReadString()+"\t"+refseq+"\t"+refseq2+"\t"+sread+"\t"+flagUnmap);
                    
//                    if(fout!=null){
                        
                     
                        if(resf.getCoverage()>=this.filter_mer_coverage&&Math.abs(pmax-pmax2)>this.minimumIndelSize&&resf.getSwitchCount()==2&&resf.getDupCount()<=this.filter_dup_count&&resf.getRefRepeatCount()<this.filter_ref_repeat){
                            int bp[] = resf.findBreakPoint();
                            String junction = read.substring(bp[0],bp[1]);
                            if(junction.length()==0)
                                junction = "-";
                            
                            String[] tA = chrposA.split(":");
                            String[] tB = chrposB.split(":");
                            int bp0 = bp[0];
                            int bp1 = bp[1];
                            String[] tpA = maxposA.split(":");
                            String[] tpB = maxposB.split(":");
                            
                            if(tpA[1].compareTo("1")==0){
                                 bp0 = sb.length()-bp0;
                            }
                            if(tpB[1].compareTo("1")==0){
                                 bp1 = sb.length()-bp1;
                            }
                            
                            
                            String junctionA = ""+tA[0]+":"+(Integer.parseInt(tA[1])+bp0);
                            String junctionB = ""+tB[0]+":"+(Integer.parseInt(tB[1])+bp1);
                            
                            
                            fout.println(""+(z+1)+"\t"+resf.getCoverage()+"\t"+resf.getSwitchCount()+"\t"+resf.getDupCount()+"\t"+maxAc+"\t"+maxBc+"\t"+maxposA+"\t"+maxposB+"\t"+chrposA+"\t"+chrposB+"\t"+r.getReadString()+"\t"+refseq+"\t"+refseq2+"\t"+sread+"\t"+resf.getResolveTag()+"\t"+bp[0]+"\t"+bp[1]+"\t"+junctionA+"\t"+junctionB+"\t"+junction+"\t"+flagUnmap);

                            
                        }
                        
                        
                        
//                    }
                    
               
                    
                    
                    }
                    
            
                    
                    
                    
                    
                    
                
                
            }
             
             
            
            
            


//            if(rx.get(max)!=null&&(int)rx.get(max)>minimum_peak_pattern1&&(int)rx.get(max2)>minimum_peak_pattern2){
//                
////              if(max2<max){
////                  long tmp = max2;
////                  max2 = max;
////                  max = tmp;
////              }  
//             
//                
//             if(debug) System.out.println(rx);
//             
//             
//             long pmax = Long.valueOf(max.split(":")[0]);
//             CombinedPos cp = this.getChromosomePos(pmax);
//             String chrA = cp.getChr();
//             int posA = cp.getPos();
//             
//             long pmax2 = Long.valueOf(max2.split(":")[0]);
//             CombinedPos cp2 = this.getChromosomePos(pmax2);
//             String chrB = cp2.getChr();
//             int posB = cp2.getPos();
//             
//             
//             String refseq, refseq2;
//             
//             if(!this.random_access){
//                refseq = this.getReferenceSequence(cp.getChr(), cp.getPos(), 150);
//                refseq2 = this.getReferenceSequence(cp2.getChr(), cp2.getPos(), 150);
//
//             }else{
//                refseq = this.getReferenceSequenceFromIndex(cp.getChr(), cp.getPos(), 150);
//                refseq2 = this.getReferenceSequenceFromIndex(cp2.getChr(), cp2.getPos(), 150);
//             }
//             
//             
//             StringBuffer sread = new StringBuffer();
//             
//             for(int x=0;x<res.size();x++){
//                 
//                 MappedResult mr = (MappedResult)res.get(x);
//                 int[] a = mr.getPos();
//                 if(a!=null){
//                     
//                     if(a.length==1){
//                         if(a[0]-x==pmax)sread.append("A");
//                         else
//                         if(a[0]-x==pmax2)sread.append("B");
//                         else
//                             sread.append(".");
//                     }
//                         
//                     
//                 }else{
//                     sread.append("_");
//                 }
//                 
//                 
//             }
//             
//             EvaluateResult resc = evaluateConfirmation(r.getReadString(),refseq,refseq2);
//             
//             if(resc.passed()){
//             
//             pout.println(""+(z+1)+"\t"+resc.getCoverage()+"\t"+resc.getSwitchCount()+"\t"+resc.getDupCount()+"\t"+rx.get(max)+"\t"+rx.get(max2)+"\t"+max+"\t"+max2+"\t"+(chrA+":"+posA)+"\t"+(chrB+":"+posB)+"\t"+r.getReadString()+"\t"+refseq+"\t"+refseq2+"\t"+sread+"\t"+resc.getResolveTag());
//             
//             
//             
//             }
//
//            }
             
         }catch(Exception e){
                StackTraceElement[] elements = e.getStackTrace();
                for(int i=0; i<elements.length; i++) {
                System.err.println(elements[i]);
                }
            }   
     //            System.out.println();
            
            
        }
        
        
        this.finishBatch(thread_id);
        
        
    }
    
    long[] encodeContig(String sb){
            
           
        
            int kmer = 16;
            int sliding = 1;
            
            int n = (sb.length()-kmer)/sliding;       
            long cmer = -1;
            long mask = 0; 
            int count =0;
            int start =0;
            
            long [] res = new long[n];

            for(int i =0;i<kmer;i++)mask=mask*4+3;


            for(int i =0;i<n;i++){

                long pos = i;
                
                
                char chx = sb.charAt(i*sliding+kmer-1);
                if(chx!='N'){
                    if(cmer==-1){
                        String s = sb.substring(i*sliding,i*sliding+kmer);
                        cmer = encodeMer(s,kmer);
                    }else{

                        int t =-1;
                        switch(chx){
                            case 'A':
                            case 'a':
                                t=0; // 00
                                break;
                            case 'T':
                            case 't': 
                                t=3; // 11
                                break;
                            case 'U':
                            case 'u': 
                                t=3; // 11 (in case of RNA T has change to U)
                                break;
                            case 'C':
                            case 'c':
                                t=1; // 01 
                                break;
                            case 'G':
                            case 'g':
                                t=2; // 10 
                                break;
                            default : 
                                t=-1;
                            break;

                        }
                        if(t>=0){

                            cmer *= 4;
                            cmer &= mask;
                            cmer += t;

                        }else{
                            cmer = -1;
                            i+=kmer;
                        }  
                    }  
                    
                    
//                    if(i%10000000==0)System.out.println("Encode - "+name+" "+(i*sliding));

//                    if(cmer>=0){ 
//                        
//                        
//                        
//                        
//                    }
                    res[i] = cmer;

                }
                
            }
                 
            return res;
                    
    }
    
    class EvaluateResult{
        
        int switch_count;
        int dup_count;
        String evaluate_tag;
        String resolve_tag = null;
        int coverage;
        String read;
        String ref1;
        String ref2;
        
        boolean reverse = false;
        
        public boolean isReversePeak(){
            return reverse;
        }
        public void setReversePeak(boolean reverse){
            this.reverse = reverse;
        }
        public boolean passed(){
            return true;
        }
        
        public int getCoverage() {
            return coverage;
        }
        
        public int getSwitchCount(){
            return switch_count;
        }
        
        public int getDupCount(){
            return dup_count;
        }
        public String getEvaluateTag(){
            return this.evaluate_tag;
        }
        public String getResolveTag(){
            return this.resolve_tag; 
        }
        
        private void setSwitchCount(int switch_count) {
            this.switch_count = switch_count;
        }

        private void setDupCount(int dup_count) {
            this.dup_count = dup_count;
        }   

        private void setResolveTag(String tag) {
            this.resolve_tag = tag;
        }
        
        private void setEvaluateTag(String tag){
            this.evaluate_tag = tag;
        }

        private void setCoverage(int coverage) {
            this.coverage = coverage;
        }

        private void setRead(String read) {
            this.read = read;
        }

        private void setRef1(String ref1) {
            this.ref1 = ref1;
        }

        private void setRef2(String ref2) {
            this.ref2 = ref2;
        }
        
        public int getRefRepeatCount(){
            int count = 0;
            
            for(int i=0;i<this.ref1.length();i++){
                char c = ref1.charAt(i);
                if(c=='a'||c=='t'||c=='c'||c=='g'){
                    count+=1;
                }
            }
            for(int i=0;i<this.ref2.length();i++){
                char c = ref2.charAt(i);
                if(c=='a'||c=='t'||c=='c'||c=='g'){
                    count+=1;
                }
            }
            return count*50/this.read.length();
            
        }
        
        private String resolve(){
            StringBuilder ntag = new StringBuilder(read);
           
            int mer = 16;
            for(int i = 0;i< this.read.length();i++){
                ntag.setCharAt(i, '_');
            }
            for(int i = 0;i< this.evaluate_tag.length();i++){
                
                char c = evaluate_tag.charAt(i);
                
                
                
                if(c=='A'||c=='B'){
                    
                    char sc = '_';
                    if(c=='A')sc = 'a';
                    if(c=='B')sc = 'b';
                    ntag.setCharAt(i, c);
                    for(int j=1;j<mer;j++){
                        ntag.setCharAt(i+j, sc);
                    }
                    
                    
                }
            }
            
            StringBuilder nntag = new StringBuilder(ntag);
            String ntagt = ntag.toString();
         
            for(int i = 0;i< ntag.length();i++){
                
                     
                if(ntagt.charAt(i)=='a'&&ntagt.charAt(i+1)=='_'){
                    
                    int j = 1;
                    while((i+j)<ntagt.length()&&ntagt.charAt(i+j)=='_'){
                        
                        if(this.read.charAt(i+j)==this.ref1.charAt(i+j)){
                            nntag.setCharAt(i+j, 'a');
                        }else{
                            nntag.setCharAt(i+j, 's');
                        }
                            
                        j++;
                    }
                    
                }
                
                if(ntagt.charAt(i)=='b'&&ntagt.charAt(i+1)=='_'){
                    
                    int j = 1;
               
                    while((i+j)<ntagt.length()&&ntagt.charAt(i+j)=='_'){
                        
                        if(this.read.charAt(i+j)==this.ref2.charAt(i+j)){
                            nntag.setCharAt(i+j, 'b');
                        }else{
                            nntag.setCharAt(i+j, 's');
                        }
//                        System.out.println(""+i+" "+j+" "+(i+j));  
                        j++;
                    }
                    
                }
                
            
            }
            
            this.setResolveTag(nntag.toString());
            
            
            return this.resolve_tag;
            
        }
        
        
        
        private int[] findBreakPoint(){
            
        
            int[] bp = new int[2];
            
            int count = 0 ;
            int last = 0 ;
            char start = ' ';
            String s = this.resolve_tag;
            
            for(int i=0;i<s.length()-1;i++){
                char c = s.charAt(i);
                char n = s.charAt(i+1);
//                System.out.println(""+start+" "+c+" "+count);
//                if(c!='_'&&n!='s'&& Character.toUpperCase(c)!=Character.toUpperCase(n)){
//                    bp[0] = i;
//                    bp[1] = i;
//                    break
//                }
                if(start==' '&&(c=='a'||c=='b'||c=='A'||c=='B')){
                    start = Character.toLowerCase(c);
                    count++;
                }else
                    if(Character.toLowerCase(c)==start)
                        count++;
                    else
                        if((start=='a'&&c=='A')||(start=='b'&&c=='B'))
                             count++;
                        else{
                            if(count>=4)
                                last = i-1;
                            count=0;
                            if(start=='a'&&c=='B'){
                                bp[0]=last;
                                bp[1]=i-1;
                                break;
                            }    
                            if(start=='b'&&c=='A'){
                                bp[0]=last;
                                bp[1]=i-1;
                                break;
                            }    
                           
                        }
                            

                
            }
 
            
            return bp;
            
        }
        
        
        
    }
    
    public EvaluateResult evaluateConfirmation(String read, String ref1, String ref2){
        int mer = 16;
        EvaluateResult r = new EvaluateResult();
        
        long[] encode_read = encodeContig(read);
        long[] encode_ref1 = encodeContig(ref1);
        
        
        if(ref2==null){
            ref2 = reverseComplement(ref1);
        }
        long[] encode_ref2 = encodeContig(ref2);
            
        
        int n = read.length()-mer;
        StringBuffer sb=new StringBuffer();
        int dup_count = 0 ;
        int switch_count = 0;
        char last = ' ';
        char lastc =  ' ';
        char next = ' ';
        int coverage = 0;
        for(int i=0;i<n;i++){
            
               if(encode_ref1[i]==encode_ref2[i]){
//               System.out.print("D "+encode_ref1[i]);
               sb.append("D");
               coverage++;
               dup_count++;
               
               }else{
               if(encode_read[i]==encode_ref1[i]){
                   sb.append("A");
                    next = 'a';
                    coverage++;
               }
               if(encode_read[i]==encode_ref2[i]){
                   sb.append("B");
                   next = 'b';
                    coverage++;
               }
              
               if(encode_read[i]!=encode_ref2[i]&&encode_read[i]!=encode_ref1[i]){
                  sb.append("_");
               }
               
               }
               
               if(next!=' '){
                   
                   if(last != next){
                       switch_count++;
                       last = next;
                       if(last!='_')
                           lastc=last;
                   }
                   
                   
               }
        
            
        }
//        System.out.println("Last C "+lastc);
        if(Character.toUpperCase(lastc)=='A'){
           r.setReversePeak(true);
        }
        
//        System.out.println();
        r.setRead(read);
        r.setRef1(ref1);
        r.setRef2(ref2);
        r.setEvaluateTag(sb.toString());
        r.setSwitchCount(switch_count);
        r.setDupCount(dup_count);
        r.setCoverage(coverage);
        return r;
    }
    
    
    
    CombineReferenceIndex index=null;
    int useMer = 1;
    boolean random_access;
    
    public int[] searchMer(int mer) throws IOException{
        
       return searchMer(mer,0);
        
    }
    
    public int[] searchMer(int mer,int thread_id) throws IOException{
        
        if(index==null){
            this.index = new CombineReferenceIndex(this);
            this.index.setRandomAccess(random_access);
        }
        return index.searchMer(mer, thread_id);
        
    }
    
    
    
    public int getMaximumDuplicatePattern(){
        return useMer;
    }
        
    public void setMaximumDuplicatePattern(int i) {
        this.useMer = i;

    }

    public void setRandomAccess(boolean b) {
           
        this.random_access = b;

    }
    
    class CombineReferenceIndex{
        
        HashMap indeies = null;
        CombineReferenceSequence ref;
        boolean loaded = false;
        boolean random_access = false;
        RandomAccessFile[] random_list ;
        
        void setRandomAccess(boolean r){
            random_access = r;
        }
        
        CombineReferenceIndex(CombineReferenceSequence ref){
            this.ref = ref;
            indeies = new HashMap();
        }
        
        public int[] searchMer(int mer) throws FileNotFoundException, IOException{
        
        return searchMer(mer,0);
        
        }
        public int[] searchMer(int mer, int thread_id) throws FileNotFoundException, IOException{
            
            int mx = ref.getMaximumDuplicatePattern();
            
            if(random_access){
                
                 if(!loaded){
                        random_list = new RandomAccessFile[ref.getMaximumDuplicatePattern()*ref.number_of_thread];
                       for(int t=0;t<ref.number_of_thread;t++)
                       for(int i=1;i<=mx;i++){
                           
                            String name = ref.filename+"."+i+".bin16.final";
                            if(i==1)name = ref.filename+".u.bin16.final";
//                             System.out.println("add "+name);
                            random_list[i-1+t*mx] = new RandomAccessFile(name,"r"); //File object
                       }
                       loaded=true;
                 }
                
                 if(loaded){
                
                 for(int i=1;i<=mx;i++){
                    
                  
//                    System.out.println(""+i+" random_list ss " +random_list);

                    RandomAccessFile fi = random_list[thread_id*mx + i-1];
//                    System.out.println(""+i+" random_list sss " +fi);
                    
                    int bin = (i+1)*4;
                    int a = 0;
                    int b = (int)(fi.length()/bin)-1;
                    int mid = 0 ;
                    long l =0;
                    long mask = mask();
                    int m=0;
//                    System.out.println(""+a+" "+b+" "+bin);
                    
                while(a<=b){
                    mid = a+ (b-a)/2;
//                     System.out.println(""+a+" "+b+" "+mid+" "+((long)mid*bin));
                    if((long)mid*bin<0)break;
                    
                    fi.seek((long)mid*bin);
//                    m = map[mid*block];
                    if(i==1)
                    m = (int)((fi.readLong()>>32)&mask);
                    else
                    m = fi.readInt();
                    
//                    System.out.println(""+a+" "+b+" "+mid+ " " +Integer.toBinaryString(m));
                    
                    
                    if(mer==m){
                        break;
                    }else
                        if(mer>m){
                            a = mid+1; 
//                            System.out.println("Gt " +m);
                        }else{
                            b = mid-1;
//                            System.out.println("Lt " +m);
                        }
                }    
              
                if(m==mer){
                    int r[] = new int[i];
                    for(int j=0;j<i;j++){
                        if(i==1){
                            fi.seek(((long)mid*bin));
                            l = fi.readLong();
                            r[0]=(int)(l&mask);
                         }    
                        else{
                                    
                            fi.seek(((long)mid*bin+4)+j*4);
                            r[j] = fi.readInt();
                        
                        }
                        
                    }
                    fi.close();
                    return r;
                }
                    
                
//                    fi.close();
                    
                 }    
                }
                
                
            return null;
            
            }else{
            if(!loaded){
                
                for(int i=1;i<=mx;i++){
                    
                    String name = ref.filename+"."+i+".bin16.final";
                    
                    if(i==1)name = ref.filename+".u.bin16.final";
                    
                    File fi = new File(name); //File object
                    DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(fi)));
                    
                    long size = fi.length();
                     System.out.println("Name "+ name+" : " + size);
                    if(i==1){
                        
                        size /= 8; 
                        long[] list = new long[(int)size];
                       
                        for(int j=0;j<size;j++){
                            
                            list[j] = is.readLong();
                            if(j%100000000==0){
                                System.out.println("loading RefIndex "+i+" : "+(j/1000000)+"M");
                            }
                            
                        }
                        
                        indeies.put(Integer.valueOf(i), list);
                        
                    }else{
                        
                        size /=4 ; 
                        int[] list = new int[(int)size];
                       
                        for(int j=0;j<size;j++){
                            
                            list[j] = is.readInt();
                            if(j%100000000==0){
                                System.out.println("loading RefIndex "+i+" : "+(j/1000000)+"M");
                            }
                            
                        }
                        
                        indeies.put(Integer.valueOf(i), list);
                    }
                    
                    
                    
                }
                
                
                
                loaded = true;
            }
            if(loaded){
                int a = 0;
                int b = 0;
                int mid=0;
                int m = 0;
                
                for(int i=1;i<=mx;i++){
                 
                
                if(i==1){
                 
                long[] map = (long[])indeies.get(i);
                
                a = 0;
                b = map.length-1;
               
                long l =0;
                long mask = mask();
                while(a<=b){
                    mid = a+ (b-a)/2;
                    l = map[mid];
                    m = (int)(l>>32);
                    if(mer==m){
                        break;
                    }else
                        if(mer>m){
                            a = mid+1; 
//                            System.out.println("Gt " +m);
                        }else{
                            b = mid-1;
//                            System.out.println("Lt " +m);
                        }
                }
                if(m==mer){
                    int r[] = {(int)(l&mask)};
                    return r;
                }
                     
                }  else{
                    
                int block = (1+i);    
                int[] map = (int[])indeies.get(i);
                
                a = 0;
                b = (map.length/block)-1;
               
                long l =0;
                long mask = mask();
                while(a<=b){
                    mid = a+ (b-a)/2;
                    m = map[mid*block];
                    if(mer==m){
                        break;
                    }else
                        if(mer>m){
                            a = mid+1; 
//                            System.out.println("Gt " +m);
                        }else{
                            b = mid-1;
//                            System.out.println("Lt " +m);
                        }
                }
                if(m==mer){
                   int r[] = new int[i];
                    for(int j=0;j<i;j++){
                        r[j] = map[mid+j];
                    }
                    
                    return r;
                }    
                    
                    
                    
                }
                   
                }
                
            }
            
            
            return null;
        }
        }
        
    }
    
    class DuplicateMer {
        int mer;  // 32 bits
        ArrayList<Integer> list; //= new ArrayList<int>();
        
        DuplicateMer(int mer){
            this.mer = mer;
            list = new ArrayList<Integer>();
        }
        
        void addPos(int pos){
            list.add(pos);
        }

        private int size() {
            return list.size();
        }

        private ArrayList<Integer> getPos() {
            return list;
        }

        private void addPos(int[] pos) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        
    }
    
    class FixedDuplicateMer implements Comparable<FixedDuplicateMer>{
        int mer;
        int pos[];
        FixedDuplicateMer(int mer, int size){
            this.mer = mer;
            this.pos = new int[size];
        }
        
        int[] getPos(){
            return pos;
        }

        @Override
        public int compareTo(FixedDuplicateMer o) {
            return this.mer-o.mer;
        }

        private void setPos(int i, int pos) {
            this.pos[i] = pos;
        }
    }
    
    
    class Worker implements Runnable{
        StringBuffer sb;
        String name;
        long list[];
        int start;
        long length;
        long offset;
        
        public Worker(String name, StringBuffer seq, long list[], int start, int length, long offset){
            sb = seq;
            this.name = name;
            this.list = list;
            this.start = start;
            this.length = length;
            this.offset = offset;
        }
        
        @Override
        public void run() {
            
            
            
            int kmer = 16;
            int sliding = 1;
            
            int n = (sb.length()-kmer)/sliding;       
            long cmer = -1;
            long mask = 0; 
            int count = 0;


            for(int i =0;i<kmer;i++)mask=mask*4+3;

            System.out.println(mask);

            for(int i =start;i<(int)(start+length)&&i<n;i++){

                long pos = i*sliding+offset;
                
                
                char chx = sb.charAt(i*sliding+kmer-1);
                if(chx!='N'){
                    if(cmer==-1){
                        String s = sb.substring(i*sliding,i*sliding+kmer);
                        cmer = encodeMer(s,kmer);
                    }else{

                        int t =-1;
                        switch(chx){
                            case 'A':
                            case 'a':
                                t=0; // 00
                                break;
                            case 'T':
                            case 't': 
                                t=3; // 11
                                break;
                            case 'U':
                            case 'u': 
                                t=3; // 11 (in case of RNA T has change to U)
                                break;
                            case 'C':
                            case 'c':
                                t=1; // 01 
                                break;
                            case 'G':
                            case 'g':
                                t=2; // 10 
                                break;
                            default : 
                                t=-1;
                            break;

                        }
                        if(t>=0){

                            cmer *= 4;
                            cmer &= mask;
                            cmer += t;

                        }else{
                            cmer = -1;
                            i+=kmer;
                        }  
                    }  
                    if(i%10000000==0)System.out.println("Encode - "+name+" "+(i*sliding));

                    if(cmer>=0){  
//                        long x = (cmer<<(64- kmer*2))|pos;
                        long x = (cmer<<32)|(pos&mask);       
                        list[i] = x;
                        count++;
                    }

                }
            }

            

        }
        
    }
    public static long mask(){
         long mask = 0; 
         for(int i =0;i<16;i++)mask=mask*4+3;
         return mask;
    }
    public void indexing() throws FileNotFoundException, IOException, InterruptedException {

        
        
        
        Enumeration<ChromosomeSequence> e = chrs.elements();
        
       // for each chromosome to each number of duplicate
       
       // allocate all bin file 
        
        
        long offset = 0 ;
        long size = 0 ;
        
        ArrayList glist = new ArrayList();
        
        
        ArrayList<String> dup_files=new ArrayList<String>();
        ArrayList<DataOutputStream> osx=new ArrayList<DataOutputStream>();
        
        for(int i=0;i<=10;i++){
            
            File f = new File(this.filename+"."+i+".bin16.part"); //File object
            DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f))); // create object for output data stream
            osx.add(i, os);
        
        }
       
        
        while(e.hasMoreElements()){
          
            ChromosomeSequence chr = e.nextElement();
//            File f = new File(chr.getFilePath()+".bin16"); //File object
//            DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f))); // create object for output data stream
//           
            

            StringBuffer sb = chr.getSequence();

            int kmer = 4;
            int sliding = 1;
            
            int n = (sb.length()-kmer)/sliding;       
            long list[] = new long[n]; // Pre - allocate Array by n
            
            int nthread = 8;
            Thread[] tlist = new Thread[nthread]; 
            
            for(int i=0;i<nthread;i++){
                int nlength = n/nthread;
                
                System.out.println(chr.getName()+" start with offset = " +offset+" in thread "+i);
                Worker w = new Worker(chr.getName(),sb,list,nlength*i,nlength,offset);
                Thread t = new Thread(w);
                t.start();
                tlist[i] = t;
            }
            
            for(int i=0;i<nthread;i++){
                 
                tlist[i].join();
            
            }
            
            
            System.out.println("Sorting..");
            Arrays.parallelSort(list);
          
            DataOutputStream osi = osx.get(1);
            long mask = mask();
            long j = -1;
            int last = -1;
            int lastpos = -1;
            DuplicateMer dup = null;
            boolean proc = false;
            
            
            
//            os.writeInt(list.length);
            for(int i=0;i<list.length;i++){
                if(i%10000000==0)System.out.println("Write  "+i);
                
                int pos = (int)(list[i]&mask);
                int mer = (int)(list[i]>>32);
                
                if(proc){
              
                
                if(mer!=last){
                    if(dup!=null){
//                        dup.addPos(lastpos);
                        // write dup mer
                        j = ((j>>32)<<32)+dup.size();
                        int s = dup.size();
                        if(s<=10){
                            
                        // write duplicate to file
                        
                        // write 0000 to 1 
                        osi.writeLong(j);
                        
                        DataOutputStream osd = osx.get(s);
                        
                        osd.writeInt(mer);
                        
                        ArrayList<Integer> l = dup.getPos();
                        for(int k=0;k<s;k++){
                             osd.writeInt(l.get(k));
                        }
                        
                        }
                        
                        
                    }else{
                        // unique write
                        osi.writeLong(j);
                    }
                    dup = null;
                }else{
                    if(dup==null){
                        dup = new DuplicateMer(mer);
                        dup.addPos(lastpos);
                    }
                    if(dup!=null){
                        dup.addPos(pos);
                    }
                }
                }
                proc = true;
//                os.writeLong(list[i]); // write list variable to file .bin
                last = mer;
                lastpos = pos;        
                j = list[i];
            
            }
//            os.close();
            
            
//            for(int i=0;i<n;i++){
//                glist.add(list[i]);
//            }
            size += list.length;
            System.out.println(chr.getName()+" "+chr.seq.length()+" "+size);
            list=null;
            System.gc();
            offset += sb.length();
         
        }     
        
        for(int i=0;i<=10;i++){
            
          
            DataOutputStream osd = osx.get(i);
            osd.close();
        
        }
       
        
        
//            Object objs[] =  glist.toArray();
//            
//            long listx[] = new long[size];
//            for(int i=0;i<objs.length;i++){
//                long l[] = (long[])objs[i];
//                
//                for(int j=0;j<l.length;j++){
//                     listx[i] = l[j];
//                }
//               
//            }
            
            
         
//        final_indexing();
            final_indexing();
            final_unique_indexing();
       
        
    }
    
    public void final_indexing() throws FileNotFoundException, IOException {
        
        
        ArrayList<String> dup_files=new ArrayList<String>();
        ArrayList<DataOutputStream> osx=new ArrayList<DataOutputStream>();
       
        for(int i=0;i<=10;i++){
            
            File f = new File(this.filename+"."+i+".bin16.part.final"); //File object
            DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f))); // create object for output data stream
            osx.add(i, os);
        
        }
        
        osx.get(0).close();
        osx.get(1).close();
        
        
        // work on 1.bin
        
        
        
        
        
        
        
        // work on other bin
        
        
        
        
        
        for(int i=1;i<=10;i++){
            DataOutputStream osc = osx.get(i);
            osc.close();
            
            
            File fi = new File(this.filename+"."+i+".bin16.part"); //File object
            DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(fi)));
           
            long size = fi.length()/((1+i)*4);
            boolean filter = false;
            boolean discard = false;
            int max_array = Integer.MAX_VALUE-8;
            int mer;
            int bid;
            
          
            
            ArrayList setlist = new ArrayList();
            
            
            int time = (int)(size/(Integer.MAX_VALUE-8))+1;
            
            System.out.println("Part "+i+" "+size+ " left "+time);
            
            if(time>1){ // large element
            
            int binc[] = new int[4];    
            for(int k=0;k<binc.length;k++)
            binc[i] = 0; 
            int binmask = 3<<28;
            
            for(int t=0;t<time;t++){
                
                int lsize = max_array;
                if(size-max_array<0)lsize = (int)(size%max_array);
                
                long list[] = new long[lsize];
                System.out.println("Part "+i+" "+size+ " left "+t+ " "+lsize+" "+max_array);
            
                for(int j=0;j<lsize;j++){
                    long l = is.readLong();
                    mer = (int)(l>>32);
                    
                    bid = 3&((int)((mer)>>30));
                    if(j%1000000==0)
                    System.out.println("Part "+i+" read "+j+" "+Integer.toBinaryString(bid));//+" "+Integer.toBinaryString(mer)+" l "+Long.toBinaryString(l));
                    binc[bid]++;
                    
                    list[j] = l;
                    
                    
                    
                }
                
                setlist.add(list);
                size -= lsize;
                
            }
            
            ArrayList outlist = new ArrayList();
            
            int a=0;
            int b=0;
            int c=0;
            int d=0;
            long al[]=null;
            long bl[]=null;
            long cl[]=null;
            long dl[]=null;
            
            
            
            for(int k=0;k<binc.length;k++){
                System.out.println("Binc "+k+" read "+binc[k]);
                if(k==0){
                    a = 0;
                    al=new long[binc[k]];
                }else
                if(k==1){
                    b = 0;
                    bl=new long[binc[k]];
                }else        
                 if(k==2){
                    c = 0;
                    cl=new long[binc[k]];
                }else
                if(k==3){
                    d = 0;
                    dl=new long[binc[k]];
                }
            }


            for(int k=0;k<setlist.size();k++){
                long plist[] = (long[])setlist.get(k);
                
                for(int j=0;j<plist.length;j++){
                    long l = plist[j];
                    mer = (int)(l>>32);
                    bid = 3&((int)((mer)>>30));
                    
                
                if(bid==0){
                    al[a++] = l;
                }else
                if(bid==1){
                    bl[b++] = l;
                }else        
                 if(bid==2){
                    cl[c++] = l;
                }else
                if(bid==3){
                    dl[d++] = l;
                }
                            
                }
                
                
            }
            
            
            outlist.add(cl);
            outlist.add(dl);
            outlist.add(al);
            outlist.add(bl);
            
            
            File fx = new File(this.filename+"."+i+".bin16.final"); //File object
            DataOutputStream osi = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fx))); // create object for output data stream
//            File fxt = new File(this.filename+"."+i+".bin16.test"); //File object
//            DataOutputStream ost = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fxt))); // create object for output data stream    
           
            
            for(int k=0;k<outlist.size();k++){
                
                long plist[] = (long[])outlist.get(k);
                
                Arrays.parallelSort(plist);
                
                
//            DataOutputStream osi = osx.get(1);
           
            long mask = mask();
            long j = -1;
            int pos = 0;
            int last = -1;
            int lastpos = -1;
            DuplicateMer dup = null;
            boolean proc = false;
            
            
            
//            os.writeInt(list.length);
            for(int n=0;n<plist.length;n++){
                if(n%10000000==0)System.out.println("Write unique part of "+i+" : "+ n);
//                ost.writeLong(plist[n]);
                pos = (int)(plist[n]&mask);
                mer = (int)(plist[n]>>32);
                
                if(proc){
                if(mer!=last){
                    if(dup!=null){
//                        dup.addPos(lastpos);
                        // write dup mer
                        int sum = 0;
                        int uni = 0;
                        int s = dup.size();
                        long luni = -1;
                        ArrayList<Integer> l = dup.getPos();

                        for(int m=0;m<s;m++){
                            
                            pos = l.get(m);
                            
                            if(pos<100)
                                sum+=pos;
                            else{
                                luni = ((j>>32)<<32)+pos;
                                sum++;
                                uni++;
                            }
                                
                        }
                        
                        j = ((j>>32)<<32)+sum;
                        
                       
                        osi.writeLong(j);
                        
                        if(uni<=10){
                        
                        if(uni>1){
                        DataOutputStream osd = osx.get(uni);
                        //&&&&&&
                        osd.writeInt(last);
                        
                        // write only uniq
                        for(int m=0;m<s;m++){
                            pos = l.get(m);
                            if(pos>100)
                            osd.writeInt(pos);
                        }
                        }else{
                           if(luni!=-1)
                           osi.writeLong(luni); 
                        }
                        }
//                        }
                        
                        
                    }else{
                        // unique write
                        osi.writeLong(j);
                    }
                    dup = null;
                }else{
                    if(dup==null){
                        dup = new DuplicateMer(mer);
                        dup.addPos(lastpos);
                    }
                    if(dup!=null){
                        dup.addPos(pos);
                    }
                }
                }
                proc = true;
//                os.writeLong(list[i]); // write list variable to file .bin
                last = mer;
                lastpos = pos;        
                j = plist[n];
            
            }
            
//            
            

            
            plist=null;
            System.gc();
            
         
        }     
         osi.close();
               
                
        }else{
            if(time==1){
                
         
            
//            osx.clear();
            
//            for(int z=0;z<=10;z++){
//            
//            File f = new File(this.filename+"."+z+".bin16.part.final2"); //File object
//            DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f,true))); // create object for output data stream
//            osx.add(z, os);
//        
//            }
            
            File fi2 = new File(this.filename+"."+i+".bin16.part.final"); //File object
            DataInputStream is2 = new DataInputStream(new BufferedInputStream(new FileInputStream(fi2)));
             
           
            //            DataOutputStream osi = osx.get(1);
            File fx = new File(this.filename+"."+i+".bin16.final"); //File object
            DataOutputStream osi = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fx))); // create object for output data stream    
            //            DataOutputStream osi = osx.get(1);
            
//            File fxt = new File(this.filename+"."+i+".bin16.test"); //File object
//            DataOutputStream ost = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fxt))); // create object for output data stream    
//            
           
            int block = (1+i);
            int is_size = (int)(fi.length()/4/block);
            int is2_size = (int)(fi2.length()/4/block);
            int to_read = is_size*block;
            int to_read2 = is2_size*block;
  
            int idx = 0;
            
          size = is_size + is2_size;
          int data[] = new int[to_read+to_read2];
            
//            size = is2_size;
//            int data[] = new int[to_read2];
            


            long idxl[] = new long[(int)size];
          
            
            System.out.println("Part "+i+" size1= "+is_size+ " size2="+is2_size);
            int v;
            
            
            // read 1
            for(int k=0;k<to_read;k++){
                v = is.readInt();
                data[k] = v;
                if(k%block==0){
                    
                    long l = ((long)v<<32)+k;
//                     System.out.println("Mer v "+Integer.toBinaryString(v)+" Mer= "+Long.toBinaryString((long)v<<32)+ " size2="+Long.toBinaryString(l));
                    
                    idxl[idx] = l;
                    idx++;
                }
                if(k%10000000==0)
                System.out.println("Read part "+i+":1 "+k);
            }
            
//            for(int z=0;z<10;z++){
//                System.out.println(Long.toBinaryString(idxl[z]));
//            }
            
            // read 2
        
            for(int k=0;k<to_read2;k++){
                v = is2.readInt();
                data[k+to_read] = v;
                if(k%block==0){
                    long l = ((long)v<<32)+(k+to_read);
                    idxl[idx] = l;
                    idx++;
                }
                if(k%10000000==0)
                System.out.println("Read part "+i+":2 "+k);
            }
            
            Arrays.parallelSort(idxl);
            
            
            long mask = mask();
            long j = -1;
            int pos = 0;
            int last = -1;
            int lastpos = -1;
            DuplicateMer dup = null;
            boolean proc = false;
                
            for(int k=0;k<idxl.length;k++){
                
                if(k%10000000==0)System.out.println("Write unique part of "+i+" : "+ k);
               
                pos = (int)(idxl[k]&mask);
                mer = (int)(idxl[k]>>32);
                
                
//                ost.writeInt(mer);
//                idx = pos;
//                for(int n=0;n<i;n++){
//                    ost.writeInt(data[idx+n+1]); 
//                }
                
                
                if(proc){
                if(mer!=last){
                    
                    
                    if(dup!=null){
//                        System.out.println("Write "+(Integer.toBinaryString(mer)+" "+Integer.toBinaryString(last))+" "+(dup.size()));
                   
                        int s = dup.size();
                        
                        ArrayList<Integer> l = dup.getPos();
                       
                       
                        if(i*l.size()<=10){
                        
//                         System.out.println("Write Dup "+(i*l.size()));
                        DataOutputStream osd = osx.get(i*l.size());
                        osd.writeInt(last);
//                        int x[] = new int[i*l.size()];
//                         System.out.println("last "+last+ " "+l.size());
                        for(int m=0;m<s;m++){
                            
                            idx = l.get(m);
//                            System.out.println("idx "+idx);
                            for(int n=0;n<i;n++){
//                                x[m*i+n] = data[idx+n];
                                osd.writeInt(data[idx+n+1]);
                                
                            }
                        }
                            
                        
                        }
//                        }
                        
                        
                    }else{
                        // unique write
//                        System.out.println("");
//                        System.out.println("Write Uniq ");

                        osi.writeInt(last);
                        idx = lastpos;
                        for(int n=0;n<i;n++){
                            osi.writeInt(data[idx+n+1]); 
                        }
                    }
                    dup = null;
                }else{
                    if(dup==null){
                        dup = new DuplicateMer(mer);
                        dup.addPos(lastpos);
                    }
                    if(dup!=null){
                        dup.addPos(pos);
                    }
                }
                }
                proc = true;
//                os.writeLong(list[i]); // write list variable to file .bin
                last = mer;
                lastpos = pos;        
                j = idxl[k];
            
            }
            
              osi.close();
            
                
            } 
      
                
        }
            
            
            
            
        }
        
        
        
    }
    //=======================================================================================================================
    //=======================================================================================================================
    //=======================================================================================================================
   
    
    
    
   public void final_unique_indexing() throws FileNotFoundException, IOException{
           File f = new File(this.filename+".1.bin16.final"); //File object
           File fo = new File(this.filename+".u.bin16.final"); //File object
            
           DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
           DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fo))); // create object for output data stream
            
            
           long size = f.length()/8;
           long mask = mask();
           for(long i=0;i<size;i++){
               
               long l = is.readLong();
               int mer = (int)(l>>32);
               int pos = (int)(l&mask);
               if(i%1000000==0)
               System.out.println("Process "+i);
               if(pos>100||pos<0){
                   os.writeLong(l);
               }
               
           }
           
           os.close();
            
            
            
            
    }
   ArrayList ref_indeies_name = null;
   ArrayList ref_indeies_seq = null;
   HashMap ref_indeies = null;
   
   public String getReferenceSequence(String chrx, int start, int length) throws FileNotFoundException, IOException{
   
        ChromosomeSequence chr = this.getChromosomeSequenceByName(chrx);
        if(start+length>chr.seq.length()){
            return null;
        }else{
            return chr.seq.substring(start, start+length);
        } 
   }

   
   
   
   public CombinedPos getChromosomePos(long pos) throws IOException{
       
       loadRefIndex();
       long ship = 0;
       for(int i=0;i<ref_indeies_seq.size();i++){
           long[] t = (long[])ref_indeies_seq.get(i);
           if(pos<t[0]){
               return new CombinedPos((String)ref_indeies_name.get(i),(int)pos);
           }else{
               pos-=t[0];
           }
       }
       return null;
   }
   
   public void loadRefIndex() throws FileNotFoundException, IOException{
       
         if(ref_indeies==null){
           
            ref_indeies = new HashMap<String,long[]>();
            ref_indeies_seq = new ArrayList<long[]>();
            ref_indeies_name = new ArrayList<String>();
            
            File index_file = new File(filename+".fai");
            String chr;
            if(index_file.exists()){
                BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(index_file)));
                while((chr=br.readLine())!=null){
               
                String[] sl = chr.split("\t");
                chr = sl[0];
//               System.out.println(chr+" "+chr.length());
                long block[] = new long[4];
                block[0] =  Long.parseLong(sl[1]);
                block[1] = Long.parseLong(sl[2]);
                block[2] =  Long.parseLong(sl[3]);
                block[3] =  Long.parseLong(sl[4]);
               
                ref_indeies.put(chr, block);
                ref_indeies_seq.add(block);
                ref_indeies_name.add(chr);
                
                
                }
                
            }else
                System.out.println("Require index file");
       }
       
       
   }
   
   public String getReferenceSequenceFromIndex(String chrx, int start, int length) throws FileNotFoundException, IOException{
       
       loadRefIndex();
       
       if(ref_indeies!=null){
           long block[] = (long[])ref_indeies.get(chrx);
           if(block!=null){
           StringBuffer sb = new StringBuffer();
           int row = (int)block[2];
           int rowi = (int)block[3];
           
           RandomAccessFile ref = new RandomAccessFile(this.filename,"r");
           int rowx = (int)(start/row);
           long s = block[1]+rowx*rowi+start%row;
           ref.seek(s);
           
           rowx = length/row+2;
           for(int i=0;i<rowx;i++){
               sb.append(ref.readLine());
           }
           
           
           
//           
//           for(int i=0;i<length;){
//               String c = Character.toString((char)ref.readByte());
//               if(s+i%rowi==0);
//               else{
//                sb.append(c);
//                i++;
//               }
//           }
           ref.close();
           return sb.substring(0, length);
           }
       }
       
       return null;
   }
   
    
   

    
}
