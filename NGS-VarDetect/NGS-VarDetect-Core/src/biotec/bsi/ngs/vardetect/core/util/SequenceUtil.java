/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.ConcatenateCut;
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
import java.math.BigInteger;
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
                
                System.out.println("No chr" +chr_file);
                
                ChromosomeSequence c = new ChromosomeSequence(ref,chr,null);
                ref.addChromosomeSequence(c);
                 
                 
                 
                File bin_file = new File(c.getFilePath()+".bin");
                File comp_bin_file = new File(c.getFilePath()+"_comp.bin");
                if(bin_file.exists()==false){
                    if(c.getSequence()==null){
                   
                    }
                     
                    EncodedSequence encoded = encodeSerialChromosomeSequenceV3(c);
                     
                     
                }
                
                if(comp_bin_file.exists()==false){
                    EncodedSequence encoded = encodeSerialChromosomeSequenceV3(c);
                }
                 
                 
            }else{
                extract_chr=true;
                System.out.println("No chr" +chr_file);

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

                        File chr_file = new File(c.getFilePath()+".fa");

                        if(!chr_file.exists()){
                            c.writeToFile("FA");
                        }                   

                        seq=null;
                        c.lazyLoad();

                        EncodedSequence encoded = encodeSerialChromosomeSequenceV3(c);

                        c.lazyLoad();

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

            File chr_file = new File(c.getFilePath()+".fa");

            if(!chr_file.exists()){
                c.writeToFile("FA");
            }     

            seq = null;

            c.lazyLoad();

            EncodedSequence encoded = encodeSerialChromosomeSequenceV3(c);

            c.lazyLoad();

            ref.addChromosomeSequence(c);
        }

    
    }catch (IOException x) {
        System.err.format("IOException: %s%n", x);
    }    
    
//    if(extract_chr){
//        Enumeration<ChromosomeSequence> e = ref.getChromosomes().elements();
//        while(e.hasMoreElements()){
//            ChromosomeSequence s = e.nextElement();
//            File chr_file = new File(s.getFilePath()+".fa");
//            if(!chr_file.exists()){
//                s.writeToFile("FA");
//            }
//        }
//    }
    
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
       
       
       File f = new File(chr.getFilePath()+".bin"); //File object
       File compF = new File(chr.getFilePath()+"_comp.bin");
       
       if(f.exists()){
           
           int count = 0 ;
            
            DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
            int size = is.readInt();
            
            long list[] = new long[size];
//            System.out.println("Totalxx bmer : "+size);

            for(int i=0;i<size;i++){
               
               
                
                list[i] = is.readLong();
                int percent = (int)(1.0*count/size); 
//                System.out.println(Long.toBinaryString(list[i]));
            
                if(percent%10==0&&percent!=0)System.out.println("Read binary Mer "+chr.getName()+" "+count);
                count ++;
            }
            
            seq.setMers(list);
//            long[] listComp = createComplimentStrand(list);
//            seq.setMersComp(listComp);
//            System.out.println("Total bmer : "+size);

            is.close();
            
            /* Compliment Part */
            if(compF.exists()){
                is = new DataInputStream(new BufferedInputStream(new FileInputStream(compF)));
                size = is.readInt();
            
                long[] complist = new long[size];
    //            System.out.println("Totalxx bmer : "+size);

                for(int i=0;i<size;i++){



                    complist[i] = is.readLong();
                    int percent = (int)(1.0*count/size); 
    //                System.out.println(Long.toBinaryString(list[i]));

                    if(percent%10==0&&percent!=0)System.out.println("Read compliment binary Mer "+chr.getName()+" "+count);
                    count ++;
                }
                seq.setMersComp(complist);
                
                complist=null;
                System.gc();
                
                is.close();
            }else{
                
                size = list.length;
                long[] complist = new long[size];
                complist = createComplimentStrand(list);
                
                list=null;
                System.gc();
                
                DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(compF)));
                os.writeInt(size);
                    for(int i=0;i<size;i++){
                        if(i%1000000==0)System.out.println("Write "+chr.getName()+" "+i);
                        os.writeLong(complist[i]); // write list variable to file .bin
                    }
                os.close();
                
                seq.setMersComp(complist);   
                
                complist=null;
                System.gc();
            }
            
//            seq.setMers(list);
   
       }else{

       DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f))); // create object for output data stream

       StringBuffer sb = chr.getSequence();
       
       
       int n = (sb.length()-kmer)/sliding;       
       long cmer = -1;
       long mask = 0; 
       int count = 0;
  
       long list[] = new long[n]; // Pre - allocate Array by n
       
              
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
            os.writeLong(list[i]); // write list variable to file .bin
        }
       os.close();
       
       seq.setMers(list);
       
       list=null;
       System.gc();
                
//       long[] listComp = createComplimentStrand(list);
//       seq.setMersComp(listComp);

    /* Compliment Part */
            if(compF.exists()){
                DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(compF)));
                int size = is.readInt();
            
                long[] complist = new long[size];
    //            System.out.println("Totalxx bmer : "+size);

                for(int i=0;i<size;i++){
                    complist[i] = is.readLong();
                    int percent = (int)(1.0*count/size); 
    //                System.out.println(Long.toBinaryString(list[i]));

                    if(percent%10==0&&percent!=0)System.out.println("Read compliment binary Mer "+chr.getName()+" "+count);
                        count ++;
                    }
                seq.setMersComp(complist);
                
                complist=null;
                System.gc();
                
                is.close();
            }else{
               
                int size = list.length;
                long[] complist = new long[size];
                complist = createComplimentStrand(list);
                
                list=null;
                System.gc();
                
                os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(compF)));
                os.writeInt(size);
                    for(int i=0;i<size;i++){
                        if(i%1000000==0)System.out.println("Write "+chr.getName()+" "+i);
                        os.writeLong(complist[i]); // write list variable to file .bin
                    }
                os.close();
                
                seq.setMersComp(complist);
                
                complist=null;
                System.gc();
            }


 
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
                  //pos<<=8;
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
    
   
    public static long encodeMer(String seq, int kmer){
       
//       cctgtagtacagtttgaagt
       long mer = 0 ;
       int t;
//       String seq = seq_input.toLowerCase();
//       if(seq.compareTo(seq_input)==0)return -1;
       for(int i=0;i < kmer && i<seq.length();i++){
            char a = seq.charAt(i);
           
            switch(a){
                case 'a':
                case 'A':
                    t=0; // 00
                    break;
                case 't': 
                case 'T':
                    t=3; // 11
                    break;
                case 'c':
                case 'C':
                    t=1; // 01 
                    break;
                case 'g':
                case 'G':
                    t=2; // 10 
                    break;
                default : 
                    t=-1;
                    return -1;
               
            }
            
            if(t>=0){

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
    
    public static EncodedSequence getEncodeSequenceV2(ChromosomeSequence chr) {

        EncodedSequence encode = null;
//        System.out.println(chr.getFilePath());
        Path fp = Paths.get(chr.getFilePath()+".bin");
        File f = fp.toFile();
        try{

            if(f.exists()){
                encode = new EncodedSequence();
                encode.readFromPath(chr.getFilePath(), "bin");
            }else{
                encode = SequenceUtil.encodeSerialChromosomeSequenceV3(chr);
    //            encode.writeToPath(chr.getFilePath(), "map");
                encode.writeToPath(chr.getFilePath(), "bin");
            }
        
        }catch(Exception e){
            System.out.println("Error YoYOYO" + e);
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
            
            System.out.println("Random cut of chr" + chrA.getChrNumber() + " from position " + iniA + " : " + cutA);
            System.out.println("Random cut of chr" + chrB.getChrNumber() + " from position " + iniB + " : " + cutB);
            System.out.println("Concatenate chromosome: " + concatenateCut);
            
            cutA = null;
            cutB = null;
            System.gc();
            
            // not finish
        //}
        
        
        return concatenateCut;
    }
    
    public static ConcatenateCut concatenateComplexChromosome(ChromosomeSequence chrA,ChromosomeSequence chrB, int cutLengthA, int cutLengthB){
        //chrA.getsequence
        int check,checkA = 1,checkB = 1;
        CharSequence cutA,cutB,dummyConcatenateCut;
        ConcatenateCut concatenateCut = new ConcatenateCut(); 
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
        concatenateCut.addBasicInfo(chrA.getName(), chrB.getName(), iniA, iniB);
        Random fusionType = new Random();
        int type = fusionType.nextInt(4);
        
        if(type == 0){
            /* strand ++ */
            dummyConcatenateCut = cutA.toString()+cutB.toString();
            concatenateCut.addSequence(dummyConcatenateCut);
            concatenateCut.addType(type);
            concatenateCut.addCutInfo(cutA, cutB);
            
            System.out.println("Random cut of chr" + chrA.getChrNumber() + " Strand (+) from position " + iniA + " : " + cutA);
            System.out.println("Random cut of chr" + chrB.getChrNumber() + " Strand (+) from position " + iniB + " : " + cutB);
            System.out.println("Concatenate chromosome: " + concatenateCut);
            
        }else if(type == 1){
            /* strand +- */
            String invCutB = SequenceUtil.inverseSequence(cutB.toString());
            String compCutB = SequenceUtil.createComplimentV2(invCutB);
            
            dummyConcatenateCut = cutA.toString()+compCutB;
            concatenateCut.addSequence(dummyConcatenateCut);
            concatenateCut.addType(type);
            concatenateCut.addCutInfo(cutA, cutB);
            
            System.out.println("Random cut of chr" + chrA.getChrNumber() + " Strand (+) from position " + iniA + " : " + cutA);
            System.out.println("Random cut of chr" + chrB.getChrNumber() + " Strand (-) from position " + (rangeB - iniB) + " : " + compCutB);
            System.out.println("Concatenate chromosome: " + concatenateCut);
            
        }else if(type == 2){
            /* strand -+ */
            String invCutA = SequenceUtil.inverseSequence(cutA.toString());
            String compCutA = SequenceUtil.createComplimentV2(invCutA);
            
            dummyConcatenateCut = compCutA + cutB.toString();
            concatenateCut.addSequence(dummyConcatenateCut);
            concatenateCut.addType(type);
            concatenateCut.addCutInfo(cutA, cutB);
            
            System.out.println("Random cut of chr" + chrA.getChrNumber() + " Strand (-) from position " + (rangeA - iniA) + " : " + compCutA);
            System.out.println("Random cut of chr" + chrB.getChrNumber() + " Strand (+) from position " + iniB + " : " + cutB);
            System.out.println("Concatenate chromosome: " + concatenateCut);
            
        }else{
            /* strand -- */
            String invCutA = SequenceUtil.inverseSequence(cutA.toString());
            String compCutA = SequenceUtil.createComplimentV2(invCutA);
            String invCutB = SequenceUtil.inverseSequence(cutB.toString());
            String compCutB = SequenceUtil.createComplimentV2(invCutB);
        
            dummyConcatenateCut = compCutA + compCutB;
            concatenateCut.addSequence(dummyConcatenateCut);
            concatenateCut.addType(type);
            concatenateCut.addCutInfo(cutA, cutB);
            
            System.out.println("Random cut of chr" + chrA.getChrNumber() + " Strand (-) from position " + (rangeA - iniA) + " : " + compCutA);
            System.out.println("Random cut of chr" + chrB.getChrNumber() + " Strand (-) from position " + (rangeB - iniB) + " : " + compCutB);
            System.out.println("Concatenate chromosome: " + concatenateCut);
            
        }
         
        cutA = null;
        cutB = null;
        System.gc();
            
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
                long roundNMatch = 0;
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
                        ///index = chrMap.get(getKeyFromValue(readMap,val)) - val;
                        ////index = chrMap.get(getKeyFromValue(readMap,val)) - roundNMatch;
                        indexB = (chrMap.get(getKeyFromValue(readMap,val))>>8) - (val>>8);
                        index = ((indexB-roundMatch)<<8)+(chrMap.get(getKeyFromValue(readMap,val))&255);
                        System.out.println("New indexB without add number of chromosome : " + indexB);
                        indexB = (indexB<<8)+(chrMap.get(getKeyFromValue(readMap,val))&255);
                        System.out.println("New indexB with add number of chromosome : " + indexB);

                        //System.out.println("Cross check with reference key : " + (getKeyFromValue(chrMap,(index+val))) );
                        System.out.println("Cross check with reference key : " + (getKeyFromValue(chrMap,(index))) );
                        //System.out.println("Position that map on reference : " + chrMap.get(getKeyFromValue(readMap,val)));

                        System.out.println("Align at : " + index );
                        //if (checkMap == 0){
                            //count++;
                        //}
                        //value[0][0] = roundMatch++;
                        //value[0][1] = read.getChrName();

                        if (roundMatch ==0 && outputMap.containsKey(index)){
                            roundMatch = 1 + (outputMap.get(index));
                            outputMap.put(index,roundMatch);
                        }
                        else{
                            roundMatch++;
                            outputMap.put(index,roundMatch);
                        }
                        
                        ///////outputMap.put(index,roundMatch);
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
    

    public static MapResult mapGenomeShotgunV4(){
        
        return null;
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

  
    public static long[] createComplimentStrand(long[] inRef){
        /*  Reconstruct whole chromosome sequence to compliment strand */
        /* inRef is contain with DNA sequence of specific chromosome which already encodeed in kmer lebgth and position */
        
        long mask = -268435456;
        long mask2 = 268435455;
        long mask36Bit = 68719476735L;
        long[] mersComp = Arrays.copyOf(inRef,inRef.length);
        
        inRef = null;
        System.gc();
        
        System.out.println("\n Create compliment strand ");
        System.out.println("******* First element before compliment : " + mersComp[0]);
       
        for(int i=0;i<mersComp.length;i++){
            long dummyMerPos = mersComp[i];
//            System.out.println("Check fullcode dummyMerPos: " + dummyMerPos);
            //long dummyMer = dummyMerPos>>28;
            long dummyMer = (dummyMerPos>>28)&mask36Bit;
            long dummyPos = dummyMerPos&mask2;

            /* Reconstruct compliment DNA sequence of whole chromosome */
//            System.out.println("Check mer sequemce befor compliment : " + dummyMer);
//            long dummyNewMer = (~dummyMer)&mask36Bit;
            
            String binaryMer = Long.toBinaryString(dummyMer);
            int kmer = binaryMer.length()/2;
//            System.out.println("Create compliment at : " + i);
//            System.out.println("Check binaryMer (before compliment) : " + binaryMer);
//            System.out.println("Check dummyPos : " + dummyPos);
            String revBin = new StringBuilder(binaryMer).reverse().toString(); //   reverse Sequence Ex 1011001 to 1001101 
//            System.out.println("Check revBin (reverse binary) : " + revBin);
            long revNum = new BigInteger(revBin,2).longValue(); //  Cast binary string to decimal number
//            System.out.println("Check revNum (decimal number of reverse binary) : " + revNum);
            long dummyNewMer = (~revNum)&mask36Bit; //  Create compliment of it
//            System.out.println("Check dummyNewMer (after compliment) : " + dummyNewMer);
            
            
            
            
//            System.out.println("Check mer sequemce after compliment : " + dummyNewMer);
//            System.out.println("Check Position before reverse : " + dummyPos);
            long dummyNewPos = (mersComp.length-1)-dummyPos; /* length-1 because assume array has 10 member ; length is 10 but maximum it index is 9 becaus index start at 0 */
            /* To get new inverse of position value, use this fomular (max index - old index) **Ex. old index is 9 so the inverse of it is (9 - 9) = 0 that's correct!! */
            // Replace to long[] 
//            System.out.println("Check Position after reverse : " + dummyNewPos);
//            System.out.println("Check max index is : " + (this.mers.length-1));
            long dummyNewMerPos = (dummyNewMer<<28)+dummyNewPos;
//            System.out.println("Check fullcode dummyNewMerPost : " + dummyNewMerPos);
            if (i%100000000==0){
                System.out.println("Create compliment at : " + i);
            }
//            System.out.println("check old code in mersComp before add " + mersComp[i]);
//            System.out.println("check code before add " + dummyNewMerPos);
            mersComp[i] = dummyNewMerPos;
//            System.out.println("check code after add " + mersComp[i]);
        }

        // re-sorted long[]
        System.out.println("******* First element after compliment unsorted  : " + mersComp[0]);
        Arrays.sort(mersComp);
        System.out.println("******* First element after compliment sorted : " + mersComp[0]);
        //return mersComp;
        return mersComp;
    }

    public static String decodeMer(long inCode, int kmer){
        
        
        String binaryCode = Long.toBinaryString(inCode);/// transform wrong
        long lengthBinary = binaryCode.length();
        long fullLen = kmer*2;
        
        if(lengthBinary<fullLen){
            for(int i =0;i<(fullLen-lengthBinary);i++){
                binaryCode = "0"+binaryCode;
            }
        }

        int lenCode = binaryCode.length();
        String sequenceBase = "";
        String dummySequenceBase = "";
//        System.out.println("this is fulLen : " + fullLen);
//        System.out.println("this is lenBin : " + lengthBinary);
//        System.out.println("this is bin code inCode : " + inCode);
//        System.out.println("this is bin code check : " + binaryCode);
        for(int i=0;i<lenCode/2;i++){
            
            int indexF = i*2;
            int indexL = indexF+2;
            String dummyBase = binaryCode.subSequence(indexF, indexL).toString();
//            System.out.println("This is base code check : " + dummyBase);
            switch(dummyBase){
                case "00":
                   dummySequenceBase = "A"; // 00     
                   break;
               case "11": 
                   dummySequenceBase = "T"; // 11
                   break;
               case "01":               
                   dummySequenceBase = "C"; // 01 
                   break;
               case "10":               
                   dummySequenceBase = "G"; // 10 
                   break;
               default : 
                   dummySequenceBase = "N";
                   //return -1;
            }
//            System.out.println("this is dummy sequence check : round = "+i+" : dummySequenceBase = " + dummySequenceBase);
            sequenceBase = sequenceBase + dummySequenceBase;    
        }
        
        return sequenceBase;
    }
    
    public static String createCompliment(String inSeq){
        /* Create compliment of short length sequence of DNA */  
        CharSequence testSeq = inSeq;
        long mask36Bit = 68719476735L;
        int kmer = testSeq.length();
        long mask = ((long)(Math.pow(2,(kmer*2)))-1) ;
        
        System.out.println("createCompliment : Check merCode = " + kmer);
        System.out.println("createCompliment : Check merCode = " + mask);
        
        long merCode = encodeMer(inSeq,kmer);
        long dummyNewMer = (~merCode)&mask;
        
        System.out.println("createCompliment : Check merCode = " + merCode);
        System.out.println("createCompliment : Check newMer = " + dummyNewMer);
        
        String outSeq = decodeMer(dummyNewMer,kmer);
        
        
        return outSeq;
    }
    
    public static String createComplimentV2(String inSeq){
        /* Create compliment of short length sequence of DNA */  
        
        int lenSeq = inSeq.length();
        String outSeq = "";
        String dummyBase;
        
        //System.out.println("this s inSeq check :\t"+inSeq);
        for (int i=0;i<lenSeq;i++){

            switch(inSeq.charAt(i)){
                case 'a':
                case 'A':
                   dummyBase = "T"; // 00     
                   break;
                case 't':   
                case 'T': 
                   dummyBase = "A"; // 11
                   break;
                case 'c':
                case 'C':               
                   dummyBase = "G"; // 01 
                   break;
                case 'g':
                case 'G':               
                   dummyBase = "C"; // 10 
                   break;
                default : 
                   dummyBase = "N";
                   //return -1;
            }
            outSeq = outSeq+dummyBase; 
        }
        //System.out.println("this s outSeq check :\t"+outSeq);
        
        
        return outSeq;
    }
    
    public static String inverseSequence(String inSeq){
        CharSequence testSeq = inSeq;
        int kmer = testSeq.length();
        String invSeq = "";
        
        for (int i = 0;i<kmer;i++){
            invSeq = invSeq + String.valueOf(testSeq.charAt((testSeq.length()-1)-i));
        }
        
        return invSeq;
    }
    
}
