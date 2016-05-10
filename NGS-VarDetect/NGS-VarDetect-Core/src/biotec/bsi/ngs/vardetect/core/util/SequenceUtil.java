/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ExonIntron;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.MapResult;
import biotec.bsi.ngs.vardetect.core.ReferenceExonIntron;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;


/**
 *
 * @author soup
 */
public class SequenceUtil {
    
//    lazy load chromosome
    
    public static ReferenceSequence getReferenceSequence(String filename) throws IOException {

    ReferenceSequence ref = new ReferenceSequence();
    ref.setFilename(filename);
    
    
//   has index then index file .index  : list of chromosome
//   read each chromosome then extract .fa to file and build bin

    System.out.println(ref.getPath());
    
    File index_file = new File(filename+".index");
    
    boolean create_index = false;
    boolean extract_chr = false;
    
    
    if(index_file.exists()){
        String chr= null;

        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(index_file)));
        while((chr=br.readLine())!=null){
            System.out.println(chr);
            File chr_file = new File(ref.getPath()+File.separator+chr+".fa");
            if(chr_file.exists()){
                 ChromosomeSequence c = new ChromosomeSequence(ref,chr,null);
                 ref.addChromosomeSequence(c);
                 
                 
                 
                 File bin_file = new File(c.getFilePath()+".bin");
                 if(bin_file.exists()==false){
                     if(c.getSequence()==null){
                   
                     }
                     EncodedSequence encoded = encodeSerialChromosomeSequenceV3(c);
                 }
                 
                 
                 
            }else{
                 extract_chr=true;
            }
        } 
    }else{
        create_index = true;
    }

    
    
    
    
    if(create_index||extract_chr){
        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(filename);
        String chr = null;
    
        StringBuffer seq = new StringBuffer();

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
        String line = null;

        while ((line = reader.readLine()) != null) {

            if(line.charAt(0)=='>'){

                if(chr!=null){

                    System.out.println("CHR : "+chr+" Size : "+seq.length());

                    ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);

                    ref.addChromosomeSequence(c);


                }
                seq = new StringBuffer();
                chr = line.substring(1,line.length());


            }else{

                seq.append(line.trim());


            }

        }
    
    if(seq.length()>0){
        
        System.out.println("CHR : "+chr+" Size : "+seq.length());
        ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);
        ref.addChromosomeSequence(c);
        
        
        
        
        
    }
    
    
    } catch (IOException x) {
        System.err.format("IOException: %s%n", x);
    }    
    
    if(extract_chr){
        Enumeration<ChromosomeSequence> e = ref.getChromosomes().elements();
        while(e.hasMoreElements()){
            ChromosomeSequence s = e.nextElement();
            File chr_file = new File(s.getFilePath()+".fa");
            if(!chr_file.exists()){
                s.writeToFile("FA");
            }
        }
    }
    
//    ========================= read all chromosome
    if(create_index){
//          File index_file = new File(ref.getPath()+".index");

        PrintStream ps = new PrintStream(new FileOutputStream(index_file));
         
        
        Enumeration<ChromosomeSequence> e2 = ref.getChromosomes().elements();
        while(e2.hasMoreElements()){
            ChromosomeSequence s = e2.nextElement();
            ps.println(s.getName());
            
        }
        
        ps.close();
    }
    
//    
    
    
    
    
    
    
    
    
    }

    
    
    return ref;
    
    }
    
    
//    public static ChromosomeSequence loadChromosomeSequence(ReferenceSequence ref, String filename);
    
    
    public static ReferenceSequence  readAndIndexReferenceSequence(String filename){
     
        
           
       
    ReferenceSequence ref = new ReferenceSequence();
    ref.setFilename(filename);
       
        
           
    Charset charset = Charset.forName("US-ASCII");
    Path path = Paths.get(filename);
    String chr = null;
//    String seq = "";
    
    StringBuffer seq = new StringBuffer();

    try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
    String line = null;
    
    while ((line = reader.readLine()) != null) {
        
        if(line.charAt(0)=='>'){
        
            if(chr!=null){
                
                System.out.println("CHR : "+chr+" Size : "+seq.length());
                
                ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);
                
                ref.addChromosomeSequence(c);
               break;
                
            }
            seq = new StringBuffer();
            chr = line.substring(1,line.length());
              
               
        }else{
            
            seq.append(line.trim());
            
            
        }
        
    }
    
    if(seq.length()>0){
        
        System.out.println("CHR : "+chr+" Size : "+seq.length());

        ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);
                
        ref.addChromosomeSequence(c);
        
    }
    
    
    
    
    System.out.println("Number of chromosomes : "+ref.getChromosomes().size());
    
    
    
    Enumeration<ChromosomeSequence> e = ref.getChromosomes().elements();
    
    
    
    while(e.hasMoreElements()){
        
        ChromosomeSequence s = e.nextElement();
        
//        EncodedSequence encoded = encodeSerialChromosomeSequenceV3(s);
        
    }
    
    
    } catch (IOException x) {
        System.err.format("IOException: %s%n", x);
    }    
           
           
     
       return ref;
        
    }
    
    
    
     public static EncodedSequence encodeSerialChromosomeSequenceV3(ChromosomeSequence chr) throws FileNotFoundException, IOException{
       
       EncodedSequence seq = new EncodedSequence();
       
       int kmer = 18;
       int sliding = 1;
       int repeat = 0;
       
       
       File f = new File(chr.getFilePath()+".bin");
       
       if(f.exists()){
           
           int count = 0 ;
            
            DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
            int size = is.readInt();
            
           long list[] = new long[size];
            System.out.println("Totalxx bmer : "+size);

            for(int i=0;i<size;i++){
               
               
                
                list[i] = is.readLong();
                int percent = (int)(1.0*count/size); 
//                System.out.println(Long.toBinaryString(list[i]));
            
                if(percent%10==0&&percent!=0)System.out.println("Read binary Mer "+chr.getName()+" "+count);
                count ++;
            }
            System.out.println("Total bmer : "+size);

            is.close();
           
            
//            seq.setMers(list);
             
           
           
           
           
       }else{
           
     
       
       DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f)));

       StringBuffer sb = chr.getSequence();
       
       
       int n = (sb.length()-kmer)/sliding;       
       long cmer = -1;
       long mask = 0; 
       int count = 0;
  
       long list[] = new long[n];
       
              
       for(int i =0;i<kmer;i++)mask=mask*4+3;
       
       System.out.println(mask);
       
       
       for(int i =0;i<n;i++){
          
          long pos = i*sliding;
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
           
          
          
           
           if(i%1000000==0)System.out.println("Encode "+chr.getName()+" "+i*sliding);
           
           if(cmer>=0){
               
              
               long x = (cmer<<(64- kmer*2))|pos;
               list[count++] = x;

           }
           
           
           }
       }
       
     
     Arrays.sort(list);

     
    os.writeInt(list.length);
    for(int i=0;i<list.length;i++){
        if(i%1000000==0)System.out.println("Write "+chr.getName()+" "+i);
        os.writeLong(list[i]);
    }
       os.close();
       
       seq.setMers(list);
 
    }
       
       
       
       
       return seq;
   }
    
    

    
    public static ReferenceSequence  readReferenceSequence(String filename){
     
       
       
    ReferenceSequence ref = new ReferenceSequence();
    ref.setFilename(filename);
       
        
           
    Charset charset = Charset.forName("US-ASCII");
    Path path = Paths.get(filename);
    String chr = null;
//    String seq = "";
    
    StringBuffer seq = new StringBuffer();

    try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
    String line = null;
    
    while ((line = reader.readLine()) != null) {
        
        if(line.charAt(0)=='>'){
        
            if(chr!=null){
                
                System.out.println("CHR : "+chr+" Size : "+seq.length());
                
                ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);
                
                ref.addChromosomeSequence(c);
                
            }
            seq = new StringBuffer();
            chr = line.substring(1,line.length());
            
        }else{
            
            seq.append(line.trim());
            
            
        }
        
    }
    
    if(seq.length()>0){
        
//        System.out.println("CHR : "+chr+" Size : "+seq.length());

        ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);
                
        ref.addChromosomeSequence(c);
        
    }
    
    
    
    
    } catch (IOException x) {
        System.err.format("IOException: %s%n", x);
    }    
           
           
     
       return ref;
   }
    

    public static void extractReferenceSequence(String filename){
     
       
       
    ReferenceSequence ref = new ReferenceSequence();
    ref.setFilename(filename);
       
        
           
    Charset charset = Charset.forName("US-ASCII");
    Path path = Paths.get(filename);
    String chr = null;
//    String seq = "";
    
    StringBuffer seq = new StringBuffer();

    try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
    String line = null;
    
    while ((line = reader.readLine()) != null) {
        
        if(line.charAt(0)=='>'){
        
            if(chr!=null){
                
                System.out.println("CHR : "+chr+" Size : "+seq.length());
                
                ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);
                
                c.writeToFile("FA");
                
                
                
            }
            seq = new StringBuffer();
            chr = line.substring(1,line.length());
            
        }else{
            
            seq.append(line.trim());
            
            
        }
        
    }
    
    if(seq.length()>0){
        
//        System.out.println("CHR : "+chr+" Size : "+seq.length());

        ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);
       
        c.writeToFile("FA");
        
    }
    
    
    
    
    } catch (IOException x) {
        System.err.format("IOException: %s%n", x);
    }    
           
           
     
   }

 
    public static EncodedSequence encodeSerialChromosomeSequence(ChromosomeSequence chr){
       
       EncodedSequence seq = new EncodedSequence();
       
       int kmer = 20;
       int sliding = 1;
       int repeat = 0;
       
       //Hashtable<Long,Long> map =new Hashtable<Long,Long>();
       //TreeMap<Long,Long> map = new TreeMap();
       Map<Long,Long> map = new HashMap();
       
       StringBuffer sb = chr.getSequence();
       String smallb = sb.toString().toLowerCase();
       
       
       
       int n = (sb.length()-kmer)/sliding;
       
       long cmer = -1;
       
       
       long mask = 0; 
       
       for(int i =0;i<kmer;i++)mask=mask*4+3;
       
       
       System.out.println(mask);
       
       
       for(int i =0;i<n;i++){
           

          String s = sb.substring(i*sliding,i*sliding+kmer);
          String s2 = smallb.substring(i*sliding, i*sliding+kmer);
          
          long pos = i*sliding;
          
          if(s.charAt(0)!='N'&&s.compareTo(s2)!=0){
              
//          System.out.println(s+" "+s2);
          if(cmer==-1){
            cmer = encodeMer(s,kmer);
          }else{
            
            char a = smallb.charAt(i*sliding+kmer-1);
            int t =-1;
            switch(a){
                case 'a':
                    t=0; // 00
                    break;
                case 't': 
                    t=3; // 11
                    break;
                case 'c':
                    t=1; // 01 
                    break;
                case 'g':
                    t=2; // 10 
                    break;
                default : 
                    t=-1;
                break;
               
            }
            if(t>=0){
                
//                String s2 = sb.substring((i-1)*sliding,(i-1)*sliding+kmer);
//                long omer = cmer;               
                cmer *= 4;
                cmer &= mask;
                cmer += t;
            
                    
//                System.out.println(""+s+" "+cmer+"\t"+s2+" "+omer+" "+t);
                
            }else{
                System.out.println(a);
                cmer = -1;
                i+=kmer;
                
                
            }
            
              
              
          }
           
          
          
           
           if(i%10000==0)System.out.println("Loop "+i*sliding+" "+repeat);
           
           if(cmer>=0){
               
               if(map.containsKey(cmer)){
                 
                  repeat ++;
                 
//                   System.out.println(s+"\t"+mer+" at "+pos+" with "+map.get(mer)); 
               }else{
                  //cmer<<=8;
                  //cmer+=chr.getChrNumber();
                  pos<<=8;
                  pos+=chr.getChrNumber();
                  map.put(cmer, pos);
               }
               
           }
           
           
           }
       }
       //System.out.println("In encodeSerialChromosome " + chrName);
       seq.setMap(map);
       
       
       System.out.println("Total mer :" +n);
       System.out.println("Total uniq mer :" +map.size());
       System.out.println("Total rep :" +repeat);
       
       
       return seq;
   }
    
    public static EncodedSequence encodeSerialChromosomeSequenceV2(ChromosomeSequence chr){
       
       EncodedSequence seq = new EncodedSequence();
       
       int kmer = 20;
       int sliding = 1;
       int repeat = 0;
       
       //Hashtable<Long,Long> map =new Hashtable<Long,Long>();
       //TreeMap<Long,Long> map = new TreeMap();
       Map<Long,Long> map = new HashMap();
       
       StringBuffer sb = chr.getSequence();
       String smallb = sb.toString().toLowerCase();
       
       
       
       int n = (sb.length()-kmer)/sliding;
       
       long cmer = -1;
       
       
       long mask = 0; 
       
       for(int i =0;i<kmer;i++)mask=mask*4+3;
       
       
       System.out.println(mask);
       
       
       for(int i =0;i<n;i++){
           

          String s = sb.substring(i*sliding,i*sliding+kmer);
          String s2 = smallb.substring(i*sliding, i*sliding+kmer);
          
          long pos = i*sliding;
          
          if(s.charAt(0)!='N'&&s.compareTo(s2)!=0){
              
//          System.out.println(s+" "+s2);
          if(cmer==-1){
            cmer = encodeMer(s,kmer);
          }else{
            
            char a = smallb.charAt(i*sliding+kmer-1);
            int t =-1;
            switch(a){
                case 'a':
                    t=0; // 00
                    break;
                case 't': 
                    t=3; // 11
                    break;
                case 'c':
                    t=1; // 01 
                    break;
                case 'g':
                    t=2; // 10 
                    break;
                default : 
                    t=-1;
                break;
               
            }
            if(t>=0){
                
//                String s2 = sb.substring((i-1)*sliding,(i-1)*sliding+kmer);
//                long omer = cmer;               
                cmer *= 4;
                cmer &= mask;
                cmer += t;
            
                    
//                System.out.println(""+s+" "+cmer+"\t"+s2+" "+omer+" "+t);
                
            }else{
                System.out.println(a);
                cmer = -1;
                i+=kmer;
                
                
            }
            
              
              
          }
           
          
          
           
           if(i%10000==0)System.out.println("Loop "+i*sliding+" "+repeat);
           
           if(cmer>=0){
               
               if(map.containsKey(cmer)){
                 
                  repeat ++;
                 
//                   System.out.println(s+"\t"+mer+" at "+pos+" with "+map.get(mer)); 
               }else{
                  
                  cmer<<=8;
                  cmer+=chr.getChrNumber();
                  
                  map.put(cmer, pos);
               }
               
           }
           
           
           }
       }
       //System.out.println("In encodeSerialChromosome " + chrName);
       seq.setMap(map);
       
       
       System.out.println("Total mer :" +n);
       System.out.println("Total uniq mer :" +map.size());
       System.out.println("Total rep :" +repeat);
       
       
       return seq;
   }
    
    public static EncodedSequence encodeSerialReadSequence(CharSequence chr){
       // Output is CharSequence
       EncodedSequence seq = new EncodedSequence();
       
       int kmer = 20;
       int sliding = 1;
       int repeat = 0;
       
       //Hashtable<Long,Long> map =new Hashtable<Long,Long>();
       //TreeMap<Long,Long> map = new TreeMap();
       Map<Long,Long> map = new HashMap();
       
       StringBuffer sb = new StringBuffer(chr);
       //StringBuffer sb = chr.getSequence();
       String smallb = sb.toString().toLowerCase();
       
       
       
       int n = (sb.length()-kmer)/sliding;
       
       long cmer = -1;
       
       
       long mask = 0; 
       
       for(int i =0;i<kmer;i++)mask=mask*4+3;
       
       
       System.out.println(mask);
       
       
       for(int i =0;i<n;i++){
           

          String s = sb.substring(i*sliding,i*sliding+kmer);
          String s2 = smallb.substring(i*sliding, i*sliding+kmer);
          
          long pos = i*sliding;
          
          if(s.charAt(0)!='N'&&s.compareTo(s2)!=0){
              
//          System.out.println(s+" "+s2);
          if(cmer==-1){
            cmer = encodeMer(s,kmer);
          }else{
            
            char a = smallb.charAt(i*sliding+kmer-1);
            int t =-1;
            switch(a){
                case 'a':
                    t=0; // 00
                    break;
                case 't': 
                    t=3; // 11
                    break;
                case 'c':
                    t=1; // 01 
                    break;
                case 'g':
                    t=2; // 10 
                    break;
                default : 
                    t=-1;
                break;
               
            }
            if(t>=0){
                
//                String s2 = sb.substring((i-1)*sliding,(i-1)*sliding+kmer);
//                long omer = cmer;               
                cmer *= 4;
                cmer &= mask;
                cmer += t;
            
                    
//                System.out.println(""+s+" "+cmer+"\t"+s2+" "+omer+" "+t);
                
            }else{
                System.out.println(a);
                cmer = -1;
                i+=kmer;
                
                
            }
            
              
              
          }
           
          
          
           
           if(i%10000==0)System.out.println("Loop "+i*sliding+" "+repeat);
           
           if(cmer>=0){
               
               if(map.containsKey(cmer)){
                 
                  repeat ++;
                 
//                   System.out.println(s+"\t"+mer+" at "+pos+" with "+map.get(mer)); 
               }else{
                  //cmer<<=8; 
                  pos<<=8;
                  map.put(cmer, pos);
               }
               
           }
           
           
           }
       }
       
       seq.setReadMap(map);
       
        
       System.out.println("Total mer :" +n);
       System.out.println("Total uniq mer :" +map.size());
       System.out.println("Total rep :" +repeat);
       
       
       return seq;
   }
  
 
    public static EncodedSequence encodeChromosomeSequence(ChromosomeSequence chr){
       
       EncodedSequence seq = new EncodedSequence();
       
       int kmer = 20;
       int sliding = 1;
       int repeat = 0;
       
       //Hashtable<Long,Long> map =new Hashtable<Long,Long>();
       //TreeMap<Long,Long> map = new TreeMap();
       Map<Long,Long> map = new HashMap();
       
       StringBuffer sb = chr.getSequence();
       
       int n = (sb.length()-kmer)/sliding;
       
       
       for(int i =0;i<n;i++){
           String s = sb.substring(i*sliding,i*sliding+kmer);
           long mer = encodeMer(s,kmer);
           long pos = i*sliding;
           
           if(i%100000==0)System.out.println("Loop "+i*sliding+" "+repeat);
           
           if(mer>=0){
               
               if(map.containsKey(mer)){
                 
                  repeat ++;
                 
//                   System.out.println(s+"\t"+mer+" at "+pos+" with "+map.get(mer)); 
               }else{
                  map.put(mer, pos);
               }
               
           }
       }
       
       seq.setMap(map);
       
       
       System.out.println("Total mer :" +n);
       System.out.println("Total uniq mer :" +map.size());
       System.out.println("Total rep :" +repeat);
       
       
       return seq;
   }
    
   
    public static long encodeMer(String seq_input, int kmer){
       
//       cctgtagtacagtttgaagt
       long mer = 0 ;
       int t;
       String seq = seq_input.toLowerCase();
       if(seq.compareTo(seq_input)==0)return -1;
       for(int i=0;i < kmer && i<seq.length();i++){
            char a = seq.charAt(i);
           
            switch(a){
                case 'a':
                    t=0; // 00
                    break;
                case 't': 
                    t=3; // 11
                    break;
                case 'c':
                    t=1; // 01 
                    break;
                case 'g':
                    t=2; // 10 
                    break;
                default : 
                    t=-1;
                    return -1;
               
            }
            
            if(t>=0){
                // ctag
                // c 1
                // ct 7
                // cta 28
                // ctag 114
                mer=mer*4+t;
                
            }
           
       }
       
       return mer;
   }

    
    public static EncodedSequence getEncodeSequence(ChromosomeSequence chr) {

        EncodedSequence encode = null;
//        System.out.println(chr.getFilePath());
        Path fp = Paths.get(chr.getFilePath()+".bmap");
        File f = fp.toFile();
        try{

        if(f.exists()){
            encode = new EncodedSequence();
            encode.readFromPath(chr.getFilePath(), "bmap");
        }else{
            encode = SequenceUtil.encodeSerialChromosomeSequence(chr);
//            encode.writeToPath(chr.getFilePath(), "map");
            encode.writeToPath(chr.getFilePath(), "bmap");
        }
        
        }catch(Exception e){
            System.out.println("Error");
        }
        //System.out.println("From Encode "+encode.getEncodeChrName());
//        return null;
        return encode;

    }
 

    public static CharSequence concatenateChromosome(ChromosomeSequence chrA,ChromosomeSequence chrB, int cutLengthA, int cutLengthB){
        //chrA.getsequence
        int check,checkA = 1,checkB = 1;
        CharSequence cutA,cutB,concatenateCut;
        CharSequence checkN = "N";
        
        int lengthA = chrA.getSequence().length();
        int lengthB = chrB.getSequence().length();
        int rangeA = lengthA - 0 ;
        int rangeB = lengthB - 0 ;
        
        Random r = new Random();
        int iniA = r.nextInt(rangeA);
        int iniB = r.nextInt(rangeB);
        cutA = chrA.getSequence().subSequence(iniA, iniA+cutLengthA);
        cutB = chrB.getSequence().subSequence(iniB, iniB+cutLengthB);
        
        
        while(checkA == 1 || checkB == 1){
            //System.out.println(rangeA);
            //System.out.println(rangeB);
            //System.out.println(lengthA);
            //System.out.println(lengthB);

            System.out.println(cutA.toString().contains("N"));
            System.out.println(cutB.toString().contains("N"));
                                
            if(cutA.toString().contains("N") || cutA.toString().equals(cutA.toString().toLowerCase())){
                iniA = r.nextInt(rangeA);
                cutA = chrA.getSequence().subSequence(iniA, iniA+cutLengthA);                
            }else checkA = 0;                        
        
            if(cutB.toString().contains("N") || cutB.toString().equals(cutB.toString().toLowerCase())){
                iniB = r.nextInt(rangeB);
                cutB = chrB.getSequence().subSequence(iniB, iniB+cutLengthB);
            }else checkB = 0;
        }
            concatenateCut = cutA.toString()+cutB.toString();
 
            System.out.println(cutA);
            System.out.println(cutB);
            System.out.println(concatenateCut);
            
            // not finish
        //}
        
        
        return concatenateCut;
    }
    
    public static Map mapGenomeShotgun(EncodedSequence refChr, InputSequence read){
        
       
        Map<Long,Long> outputMap = new HashMap();
        Map<Long,Long> chrMap = refChr.getEncodeMap();
        Long[][] value = new Long[1][2];
        String name;
                
        //Long chrnum = refChr.getEncodeChrName().split("chr");
        Vector readVector = read.getInputSequence();
        Enumeration<ShortgunSequence> e = readVector.elements();
        while(e.hasMoreElements()){
            ShortgunSequence ss = e.nextElement();
          
            EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(ss.getSequence());
            
            SortedSet<Long> vals = new TreeSet(encodeSim.getEncodeMap().values()) ;
            Map<Long,Long> readMap = encodeSim.getEncodeMap();
         
            int count = 1;
            long index = 0;
            long indexB = 0;
            long roundMatch = 1;
            long roundNMatch = 1;
            //Obiect readKey = 0;
            for (Long val : vals){
                //readKey = getKeyFromValue(readMap,val);
                //System.out.println(key);
                System.out.println("Number of mapping iteration : " + count++);
    //            System.out.println("the value is " + val + " Key is " + getKeyFromValue(read.getEncodeMap(),val));
                if (chrMap.containsKey(getKeyFromValue(readMap,val))){

                    //System.out.println("Key " + getKeyFromValue(read.getEncodeMap(),val) +": Match at position " + chr.getEncodeMap().get(getKeyFromValue(read.getEncodeMap(),val)));
                    System.out.println("val : " + val);
                            
                    System.out.println("Key " + getKeyFromValue(readMap,val) +": Match at position " + chrMap.get(getKeyFromValue(readMap,val)));
                    index = chrMap.get(getKeyFromValue(readMap,val)) - val;
                    
                    indexB = (chrMap.get(getKeyFromValue(readMap,val))>>8) - (val>>8);
                    System.out.println("New indexB without add number of chromosome : " + indexB);
                    indexB = (indexB<<8)+(chrMap.get(getKeyFromValue(readMap,val))&255);
                    System.out.println("New indexB with add number of chromosome : " + indexB);
                     
                    System.out.println("Cross check with reference key : " + (getKeyFromValue(chrMap,(index+val))) );
                    //System.out.println("Position that map on reference : " + chrMap.get(getKeyFromValue(readMap,val)));
                    
                    System.out.println("Align at : " + index );
                    //if (checkMap == 0){
                        //count++;
                    //}
                    //value[0][0] = roundMatch++;
                    //value[0][1] = read.getChrName();
                    outputMap.put(index,roundMatch++);
                     

                }else{
                    System.out.println(" Key " + getKeyFromValue(readMap,val) + ": Not match");
                    //index = val;
                    roundMatch = 1;  // reset round
                    //System.out.println("");
                }   

            }

            System.out.println("Key : "+index+ "Has : "+outputMap.get(index));
        }
        
        
        
        return outputMap;
    }
    
    
    public static MapResult mapGenomeShotgunV2(EncodedSequence refChr, InputSequence read){
        
        MapResult result = new MapResult();
        Map<Long,Long> outputMap = new HashMap();
        Map<Long,Long> chrMap = refChr.getEncodeMap();
        Long[][] value = new Long[1][2];
        String name;
                
        //Long chrnum = refChr.getEncodeChrName().split("chr");
        Vector readVector = read.getInputSequence();
        Enumeration<ShortgunSequence> e = readVector.elements();
        while(e.hasMoreElements()){
            ShortgunSequence ss = e.nextElement();
          
            EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(ss.getSequence());
            
            SortedSet<Long> vals = new TreeSet(encodeSim.getEncodeMap().values()) ;
            Map<Long,Long> readMap = encodeSim.getEncodeMap();
         
            int count = 1;
            long index = 0;
            long indexB = 0;
            long roundMatch = 1;
            long roundNMatch = 1;
            //Obiect readKey = 0;
            for (Long val : vals){
                //readKey = getKeyFromValue(readMap,val);
                //System.out.println(key);
                System.out.println("Number of mapping iteration : " + count++);
    //            System.out.println("the value is " + val + " Key is " + getKeyFromValue(read.getEncodeMap(),val));
                if (chrMap.containsKey(getKeyFromValue(readMap,val))){

                    //System.out.println("Key " + getKeyFromValue(read.getEncodeMap(),val) +": Match at position " + chr.getEncodeMap().get(getKeyFromValue(read.getEncodeMap(),val)));
                    System.out.println("val : " + val);
                            
                    System.out.println("Key " + getKeyFromValue(readMap,val) +": Match at position " + chrMap.get(getKeyFromValue(readMap,val)));
                    index = chrMap.get(getKeyFromValue(readMap,val)) - val;
                    
                    indexB = (chrMap.get(getKeyFromValue(readMap,val))>>8) - (val>>8);
                    System.out.println("New indexB without add number of chromosome : " + indexB);
                    indexB = (indexB<<8)+(chrMap.get(getKeyFromValue(readMap,val))&255);
                    System.out.println("New indexB with add number of chromosome : " + indexB);
                     
                    System.out.println("Cross check with reference key : " + (getKeyFromValue(chrMap,(index+val))) );
                    //System.out.println("Position that map on reference : " + chrMap.get(getKeyFromValue(readMap,val)));
                    
                    System.out.println("Align at : " + index );
                    //if (checkMap == 0){
                        //count++;
                    //}
                    //value[0][0] = roundMatch++;
                    //value[0][1] = read.getChrName();
                    /////outputMap.put(index,roundMatch++);
                    result.addResultMap(index,roundMatch++);
                     

                }else{
                    System.out.println(" Key " + getKeyFromValue(readMap,val) + ": Not match");
                    //index = val;
                    roundMatch = 1;  // reset round
                    //System.out.println("");
                }   
                
            }

            System.out.println("Key : "+index+ "Has : "+result.getResultMap().get(index));
        }
        
        
        
        return result;
    }
    
 public static MapResult mapGenomeShotgunV3(ReferenceSequence refGene, InputSequence read,int startElement, int stopElement) throws IOException{
        
        // Process each read separately and add result to array list
        // chr 21,and 22 is at element 12,13 of whole genome reference
        MapResult result = new MapResult();
        Vector readVector = read.getInputSequence();
        
        for(int numElement = startElement;numElement<=stopElement;numElement++){
            ChromosomeSequence chr = refGene.getChromosomes().elementAt(numElement);
            System.out.println("Name of Pick chromosome : " + chr.getName());

            EncodedSequence encode = SequenceUtil.getEncodeSequence(chr);
            System.out.println("Size of encode referecce : " + encode.getEncodeMap().size());

            //Long chrnum = refChr.getEncodeChrName().split("chr");
            
            
            //Vector readVector = read.getInputSequence();
            Enumeration<ShortgunSequence> e = readVector.elements();
            while(e.hasMoreElements()){

                Map<Long,Long> outputMap = new HashMap();
                Map<Long,Long> chrMap = encode.getEncodeMap();
                Long[][] value = new Long[1][2];
                String name;




                ShortgunSequence ss = e.nextElement();

                EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(ss.getSequence());

                SortedSet<Long> vals = new TreeSet(encodeSim.getEncodeMap().values()) ;
                Map<Long,Long> readMap = encodeSim.getEncodeMap();

                int count = 1;
                long index = 0;
                long indexB = 0;
                long roundMatch = 0;
                long roundNMatch = 1;
                //Obiect readKey = 0;

                for (Long val : vals){
                    //readKey = getKeyFromValue(readMap,val);
                    //System.out.println(key);
                    System.out.println("Number of mapping iteration : " + count++);
        //            System.out.println("the value is " + val + " Key is " + getKeyFromValue(read.getEncodeMap(),val));
                    if (chrMap.containsKey(getKeyFromValue(readMap,val))){

                        //System.out.println("Key " + getKeyFromValue(read.getEncodeMap(),val) +": Match at position " + chr.getEncodeMap().get(getKeyFromValue(read.getEncodeMap(),val)));
                        System.out.println("val : " + val);

                        System.out.println("Key " + getKeyFromValue(readMap,val) +": Match at position " + chrMap.get(getKeyFromValue(readMap,val)));
                        index = chrMap.get(getKeyFromValue(readMap,val)) - val;

                        indexB = (chrMap.get(getKeyFromValue(readMap,val))>>8) - (val>>8);
                        System.out.println("New indexB without add number of chromosome : " + indexB);
                        indexB = (indexB<<8)+(chrMap.get(getKeyFromValue(readMap,val))&255);
                        System.out.println("New indexB with add number of chromosome : " + indexB);

                        System.out.println("Cross check with reference key : " + (getKeyFromValue(chrMap,(index+val))) );
                        //System.out.println("Position that map on reference : " + chrMap.get(getKeyFromValue(readMap,val)));

                        System.out.println("Align at : " + index );
                        //if (checkMap == 0){
                            //count++;
                        //}
                        //value[0][0] = roundMatch++;
                        //value[0][1] = read.getChrName();

                        if (roundMatch ==0 && outputMap.containsKey(index)){
                            roundMatch = roundMatch + (outputMap.get(index));
                        }
                        else{
                            roundMatch++;
                        }
                        outputMap.put(index,roundMatch);
                        /////result.addResultMap(index,roundMatch++);


                    }else{
                        System.out.println(" Key " + getKeyFromValue(readMap,val) + ": Not match");
                        //index = val;
                        roundMatch = 0;  // reset round
                        //System.out.println("");
                    }   

                }

                //result.addResultArray(outputMap);

                System.out.println("Key : "+index+ "Has : " + outputMap.get(index));
                result.addResult(outputMap,ss.getReadName());

            }
        }
        
        /* Save result to file */
        
        //result.writeToPath(chr.getFilePath(),"map");
        result.writeToPath(refGene.getPath(), "map");
        return result;
    }    
    
    
    /*public static MapResult mapping(ReferenceSequence refGene, InputSequence read){
        
        MapResult result = new MapResult();
        Map outputMap;
        Enumeration<ShortgunSequence> e = read.seqs.elements();
        while(e.hasMoreElements()){
            
            ShortgunSequence ss = e.nextElement();
            outputMap = SequenceUtil.mapGenomeShotgunV3(refGene, ss);
            result.addResultArray(outputMap);
        }
        
        return result;
    }*/
    
    public static Map mapGenome(EncodedSequence chr, EncodedSequence read){
        
        
        //TreeMap<Long,Long> ref = new TreeMap(chr.getEncodeMap());
        //TreeMap<Long,Long> test = new TreeMap(read.getEncodeMap());
        SortedSet<Long> vals = new TreeSet(read.getEncodeMap().values()) ;
        Map<Long,Long> outputMap = new HashMap();
        Map<Long,Long> chrMap = chr.getEncodeMap();
        Map<Long,Long> readMap = read.getEncodeMap();
        MapResult result = new MapResult();
        
        int count = 1;
        long index = 0;
        long roundMatch = 1;
        long roundNMatch = 1;
        long dummyPos;
        long dummyIndex;
        long dummyChrName;
        Object dummykey;
        //Obiect readKey = 0;
        for (Long val : vals){
            //readKey = getKeyFromValue(readMap,val);
            //System.out.println(key);
            System.out.println("Number of mapping iteration : " + count++);
//            System.out.println("the value is " + val + " Key is " + getKeyFromValue(read.getEncodeMap(),val));
            dummykey = getKeyFromValue(readMap,val);
            
            //if (chrMap.containsKey(getKeyFromValue(readMap,val))){
            if (chrMap.containsKey(dummykey)){    
                //System.out.println("Key " + getKeyFromValue(read.getEncodeMap(),val) +": Match at position " + chr.getEncodeMap().get(getKeyFromValue(read.getEncodeMap(),val)));
                //System.out.println("Key " + getKeyFromValue(chrMap,val) +": Match at position " + chrMap.get(getKeyFromValue(readMap,val)));
                System.out.println("Key " + dummykey +": Match at position " + chrMap.get(dummykey));
                //index = chrMap.get(getKeyFromValue(readMap,val)) - val;
                
                /* Do more OOP for convert position and get chr number
                    Now I put chromosome number in position */
                
                
                
                dummyPos = chrMap.get(dummykey)>>8;
                dummyChrName = chrMap.get(dummykey)&255;
                index = dummyPos - (val>>=8);
                index = (index<<=8)+dummyChrName;
                System.out.println("Align at : " + index );
                //if (checkMap == 0){
                    //count++;
                //}
                outputMap.put(index,roundMatch++);
                
                
            }else{
                System.out.println(" Key " + getKeyFromValue(readMap,val) + ": Not match");
                //index = val;
                roundMatch = 1;  // reset round
                //System.out.println("");
            }   
            
        }
        
        System.out.println("Key : "+index+ "Has : "+outputMap.get(index)); // Print the number of align at index variable(Key) 
                                                                            // outputMap is hashmap the store the align index (Key) and number of base that align on that index.
        //TreeMap<Long> ref = new TreeMap(chr.);
        //System.out.println(test);
        //Enumeration e = read.getEncodeMap().keys();
        
        //Iterator e = test.;
        
        /*for (Map.Entry<Long,Long> entry : test.entrySet()){
            Long key = entry.getKey();
            
            if (ref.containsKey(key)){
                System.out.println("Key " + key +": Match at position " + ref.get(key));
            }else{
                System.out.println("Key " + key + ": Not match");
            }   
        }*/
            
            
            
       
        /*while(e.hasMoreElements()){
            
        
        //for (int i=0 ;i<read.getEncodeMap().size();i++){
            
            //System.out.println(e.nextElement());
            //System.out.println(chr.getEncodeMap().containsKey(e.nextElement()));
            if (chr.getEncodeMap().containsKey(e.nextElement())){
                System.out.println("Key " + e.nextElement()+": Match at position " +chr.getEncodeMap().get(e.nextElement()));
            }else{
                System.out.println("Key " + e.nextElement()+ ": Not match");
            }
        }*/
        return outputMap;
    }

    
    private static Object getKeyFromValue(Map hm, Object value) {
        for (Object o : hm.keySet()) {
          if (hm.get(o).equals(value)) {
            return o;
          }
        }
        return null;    
    }
    
    
    public static ReferenceExonIntron readExonIntron(String filename){
        
        ReferenceExonIntron ref = new ReferenceExonIntron();
        //ref.setFilename(filename);
        //Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(filename);
        String chr = null;
    //    String seq = "";
        StringBuffer seq = new StringBuffer();
        
        try (BufferedReader reader = new BufferedReader(new FileReader(filename));) {
        String line = null;
        String setA[] = null;
        String setB[] = null;
        String chrName = null;
        String geneName = null;
        long startPos = 0;
        long stopPos = 0;
        int direction = 0;
        
        while ((line = reader.readLine()) != null) {
            
            int count = 0;
            setA = line.split("\\s+");
            setB = line.split(";");
            //System.out.println(line);
            //System.out.println(setA.length);
            
            for (String part : setA) {
                count++;
                if (count == 1){
                    chrName = part;
                }
                else if(count == 4){
                    startPos = Long.parseLong(part);
                    //System.out.println(startPos);
                }
                else if (count == 5){
                    stopPos = Long.parseLong(part);
                }
                else if (count ==6){
                    direction = Integer.parseInt(part);
                }
                
                //System.out.println(part);
            }
            
            count = 0;
            for (String part : setB) {
                count++;
                if (count ==3){
                    geneName = part;
                }
                //System.out.println(part);
            }
            
            System.out.println("Get :" + chrName+ "  "  + geneName+ "  " + startPos+ "  " + stopPos+ "  " + direction);
            
            ExonIntron data = new ExonIntron(chrName,geneName,startPos,stopPos,direction);
            ref.addData(data);
            
//System.out.println(line.split(";"));
            /*if(line.charAt(0)=='>'){

                if(chr!=null){

                    System.out.println("CHR : "+chr+" Size : "+seq.length());

                    ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);

                    ref.addChromosomeSequence(c);

                }
                seq = new StringBuffer();
                chr = line.substring(1,line.length());

            }else{

                seq.append(line.trim());


            }*/

        }

        /*if(seq.length()>0){

    //        System.out.println("CHR : "+chr+" Size : "+seq.length());

            ChromosomeSequence c = new ChromosomeSequence(ref,chr,seq);

            ref.addChromosomeSequence(c);

        }*/




        } catch (IOException x) {
            System.err.format("IOException: %s%n", x);
        }    


        return ref;
    }
    
    
    public static ReferenceExonIntron randomExonIntron(ReferenceExonIntron ref){
        ReferenceExonIntron output = new ReferenceExonIntron();
        Vector<ExonIntron> data = ref.getdata();
        Random rand = new Random();
        //Vector<ExonIntron> output;
        int r = rand.nextInt(data.size());
        System.out.println("Vector size" + data.size());
        
        //int r = 500;
        System.out.print("Pick" + data.elementAt(r).getdirection());
        int iniguess = Math.abs(data.elementAt(r).getdirection());
        int iniguesscheck = Math.abs(data.elementAt(r).getdirection());
        System.out.println("iniGuesscheck" + iniguess);
        int mainPoint = 0;
        int axonNum = 0;
        int rguess = r;
        
        while(true){
            rguess++;
            System.out.println("num of r : "+rguess);
            int newguess = Math.abs(data.elementAt(rguess).getdirection());
            System.out.println("num of guess : "+newguess);
            if(newguess<iniguess||newguess==iniguess){
                mainPoint = rguess-1;
                axonNum = iniguess;
                break;
            }
            else{
                iniguess = newguess;
            }
        }
        System.out.println("mainPoint : "+ mainPoint);
        System.out.println("axonNum : "+ axonNum);
        
        String name = data.elementAt(mainPoint).getGeneName();
        //output = new Vector<ExonIntron>();
        while(axonNum>0){
            
            
            if (data.elementAt(mainPoint).getGeneName().equalsIgnoreCase(name)){
                // check in case that direction of exon is more than 2 but have just one exon EX. around element 1000
                System.out.println("At Point : " + mainPoint);
                System.out.println("Final Pick " + data.elementAt(mainPoint).getGeneName()+" Direction" + data.elementAt(mainPoint).getdirection());
                //output.add(data.elementAt(mainPoint));
                output.addData(data.elementAt(mainPoint));
                axonNum--;
                mainPoint--;
            }else break;
            
            //axonNum--;
            //mainPoint--;
            //r++;
        }
        
        //System.out.println("New vector size : "+output.getdata().size());
        output.reverseorder();
        
        //System.out.println("New vector size : "+output.getdata().elementAt(0).getdirection());
        //System.out.println("New vector size : "+output.getdata().elementAt(1).getdirection());
        
        return output;
    }

    public static void createHistrogram(MapResult inResult, String readName){
        ArrayList indexArray = new ArrayList();
        
        indexArray = indexOfAll(readName,inResult.getReadName());
        
        System.out.println("\tHistogram of "+readName);
        
        for (int i = 0;i<indexArray.size();i++){
            int index = Integer.valueOf(indexArray.get(i).toString());
            long Position = Long.valueOf(inResult.getAlignPosition().get(index).toString());
            String readN = inResult.getReadName().get(index).toString();
            long match = Long.valueOf(inResult.getNumMatch().get(index).toString());
            long chrName = Long.valueOf(inResult.getchrNumber().get(index).toString().toString());
            String output = "";
            
            for (int j=0;j<match;j++){
                output +="*";
            }
            
            System.out.println(readN+"\t|\tchr"+chrName+"\t|\t"+Position+"\t| "+output);
            
        }     
    }
    


    public static ArrayList indexOfAll(String obj, ArrayList list){
        ArrayList indexList = new ArrayList();
        for (int i = 0; i < list.size(); i++)
            if(obj.equals(list.get(i)))
                indexList.add(i);
        return indexList;
    }

  


}