/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author soup
 */
public class EncodedSequence {

    //Hashtable<Long,Long> map;
    //TreeMap<Long,Long> map;
    long mers[];            // Store reference sequence of specific chr for mapping propose
    
    Map<Long,Long> map;
    String name;
    
    public long mask = -268435456;
    public long mask2 = 268435455;
    public long mask36Bit = 68719476735L;
    
    
    void EncodedSequence(){
    }
    
    public long[] align2(long mer){
        
        int index = align(mer, 0, mers.length-1);
        
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
                    j[i-start] = mers[i]&mask2;
                    System.out.println();
                    System.out.println("Check mers the value should be 64 bit : " + mers[i] );
                    System.out.println("Check j the value should be 28 bit : " + j[i-start]);
                    System.out.println();
                }


                return j;
            }
            else
                return null;
//            System.out.println("start : "+start+" stop : "+stop+" length :"+(stop-start));
            
            
            
        }
        
        createComplimentStrand();
        return null;
    }
    
    public long align(long mer){
        
        int index = align(mer, 0, mers.length-1);
        if(index>0)return mers[index]&mask2;
        return -1;
    }
    
    public int align(long mer, int left, int right){
        
        int mid = (left+right)/2;
        long i = mers[mid]&mask;
        
        if(left>right)return -1;
        else
            if(i<mer){
                return align(mer, mid+1,right);
            }else
                if(i>mer){
                    return align(mer, left,mid-1);
                }else
                    if(i==mer){
//                        long j = mers[mid];
                        return mid;
                    }
        
        return -1;
    }
    
    
    
    
    
    
    
    
    
    public void setMers(long mers[]){
        this.mers = mers;     
    }
    
    public long [] getMers(){
        return this.mers;
    }
    
    
    public void setMap(Map<Long, Long> map) {
        this.map = map;
        //this.name = chrName;
    }
    
    public void setReadMap(Map<Long,Long> map){
        this.map = map;
    }
    
    public Map getEncodeMap(){
        return this.map;
    }
    
    public String getEncodeChrName(){
        return this.name;
    }
    
    public void readFromPath(String file_path, String fa) throws FileNotFoundException, IOException {
        
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
    
    
    public void writeToPath(String path, String fa) throws FileNotFoundException, IOException {

       
       
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

    }
    
    public void createComplimentStrand(){
        
        System.out.println(" Create compliment strand ");
        for(int i=0;i<this.mers.length;i++){
            long dummyMerPos = this.mers[i];
            System.out.println("Check fullcode dummyMerPos: " + dummyMerPos);
            long dummyMer = dummyMerPos>>28;
            long dummyPos = dummyMerPos&mask2;
            
            // Reconstruct (compliment DNA sequence)
            System.out.println("Check mer sequemce befor compliment : " + dummyMer);
            long dummyNewMer = (~dummyMer)&mask36Bit;
            System.out.println("Check mer sequemce after compliment : " + dummyNewMer);
            System.out.println("Check Position before reverse : " + dummyPos);
            long dummyNewPos = (this.mers.length-1)-dummyPos; // length-1 because assume array has 10 member ; length is 10 but maximum it index is 9 becaus index start at 0 
            // To get new inverse of position value, use this fomular (max index - old index) **Ex. old index is 9 so the inverse of it is (9 - 9) = 0 that's correct!! **
            // Replace to long[] 
            System.out.println("Check Position after reverse : " + dummyNewPos);
            System.out.println("Check max index is : " + (this.mers.length-1));
            long dummyNewMerPos = (dummyNewMer<<28)+dummyNewPos;
            System.out.println("Check fullcode dummyNewMerPost : " + dummyNewMerPos);
        }
        
        // re-sorted long[]
        Arrays.sort(this.mers);
        
    }
               
               

}
