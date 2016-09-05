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
    long[] mers;            // Store reference sequence of specific chr for mapping propose
    long[] mersComp;
    
    Map<Long,Long> map;
    String name;
    
    public long mask = -268435456; // Do & operation to get mer  
    public long mask2 = 268435455; // Do & operation to get position
    public long mask36Bit = 68719476735L;
    
    
    public EncodedSequence(){
        this.mers = null;            // Store reference sequence of specific chr for mapping propose
        this.mersComp = null;
    }
    
    public long addStrandNotation(long in, int strand){
        
        long pos_strand = (strand<<28)+in;      // add 1 bit to indicate strand notation 0(-) or 1(+) infront of the 28 position bit
        
        return pos_strand;
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
                    j[i-start] = addStrandNotation(mers[i]&mask2,strand);
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
    
    public long align(long mer){        // Curectly unuse
        
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
               
               

}
