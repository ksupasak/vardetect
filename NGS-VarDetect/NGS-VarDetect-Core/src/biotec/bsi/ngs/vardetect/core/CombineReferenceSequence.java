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
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author soup
 */
public class CombineReferenceSequence extends ReferenceSequence{
    
    
    String output = null;
    PrintStream pout =null;
    
    long skip_read = 0;
    long total_read = Integer.MAX_VALUE;
    
    int number_of_thread = 1;
    
    int minimum_peak_pattern1 = 10;
    int minimum_peak_pattern2 = 5;
    
    
    
    
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
        
        boolean finished=false;
        
        pool = new ThreadPool(number_of_thread);
        long time = this.total_read/number_of_thread;
        
        for (int i = 0; i < number_of_thread ; i++) {
            SVProfiler task = new SVProfiler(this,seq_file ,time,time*i+skip_read, i);
            queue.put(i, task);
            pool.execute(task);
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
    
    
    
    
    public  void profileSV(String filename,long countr,long skip, int thread_id) throws IOException{
        
        
        
        
        File bamFile = new File(filename); 
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamReader samReader = samReaderFactory.open(bamFile);

        
        
        SAMRecordIterator rr = samReader.iterator();
        int last_percent = 0 ;
        for(long z=0;z<countr+skip-1;z++){
            
            int percent = (int)((z-skip)*100/countr);
            
            if(percent>0&&percent!=last_percent){        
                System.out.println("Process "+thread_id+ " "+percent+"%");
                last_percent = percent;
            }
            
            SAMRecord r = rr.next();
            if(z<skip-1)continue;
            
            String sb = r.getReadString();
            
              
            int kmer = 16;
            int sliding = 1;
            
            int n = (sb.length()-kmer)/sliding;       
            long cmer = -1;
            long mask = 0; 
            int count =0;
            int start =0;

            HashMap rx = new HashMap();
            
            
            long max = 0 ;
            long max2 = 0 ;
            
            boolean debug=true;
            
            
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

                    if(cmer>=0){  
//                        long x = (cmer<<(64- kmer*2))|pos;
//                        long x = (cmer<<32)|pos;        // left Shift 28 bit is mean we reserve 28 bit for position (we fixed it) 
//                        list[i] = x;
                        int[] result = this.searchMer((int)cmer,thread_id);
                        
                        
                        
                        if(result!=null){
//                            System.out.print(""+i+ " "+result.length);
                            
                            for(int j=0;j<result.length;j++){
                            
                                pos = (((long)result[j])&mask);
                                
                                if(debug)System.out.print(" "+pos);
                                
                                
                                long xpos = pos - i;
                                if(rx.containsKey(xpos)){
                                    int v = (int)rx.get(xpos)+1;
                                    rx.put(xpos,v);
                                    if(rx.get(max)==null)max=xpos;
                                    if(v>(int)rx.get(max)){
                                        max = xpos;
                                    }
                                           
                                    
                                    
                                    
                                    
                                }else
                                   rx.put(xpos, 1);
                                
                            }
                            
                             if(debug)System.out.println();
                            
//                            for(int j=0;j<result.length;j++)  System.out.print(" "+result[j]);
                            
//                            System.out.println();
                            
//                            System.out.print(result.length);
                          
                            
                        }else{
                             if(debug)System.out.println("-");
                        }
                        
                        

                        count++;
                    }
                    
                

                }
                
              
                 
            }

//            rx.entrySet()
            
            Iterator it = rx.entrySet().iterator();
//            m
            while (it.hasNext()) {
                    
                    Map.Entry pair = (Map.Entry)it.next();
                    long k = (long)pair.getKey();
                    int v = (int)pair.getValue();
                    if(k!=max){
                        if(rx.get(max2)==null)max2=k;
                        else{
                            if(v>(int)rx.get(max2)){
                                max2=k;
                            }
                        }
                        
                    }
//                    it.remove(); // avoids a ConcurrentModificationException
            }
            
            if(rx.get(max)!=null&&(int)rx.get(max)>minimum_peak_pattern1&&(int)rx.get(max2)>minimum_peak_pattern2){
                
                
             if(debug) System.out.println(rx);
             
             pout.println(""+(z+1)+"\t"+rx.get(max)+"\t"+rx.get(max2)+"\t"+max+"\t"+max2+"\t"+r.getReadString());
                
                   

            }
             
            
     //            System.out.println();
            
            
        }
        
        
        this.finishBatch(thread_id);
        
        
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
                        long x = (cmer<<32)|pos;        // left Shift 28 bit is mean we reserve 28 bit for position (we fixed it) 
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
        int size = 0 ;
        
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
            
            int nthread = 2;
            Thread[] tlist = new Thread[nthread]; 
            
            for(int i=0;i<nthread;i++){
                int nlength = n/nthread;
                
                
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
               if(pos>100){
                   os.writeLong(l);
               }
               
           }
           
           os.close();
            
            
            
            
    }

    
}
