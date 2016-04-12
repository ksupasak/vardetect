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
import biotec.bsi.ngs.vardetect.core.ReferenceExonIntron;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
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
                  map.put(cmer, pos);
               }
               
           }
           
           
           }
       }
       
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
                  map.put(cmer, pos);
               }
               
           }
           
           
           }
       }
       
       seq.setMap(map);
       
        
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
    
    
    public static Map mapGenome(EncodedSequence chr, EncodedSequence read){
        
        
        //TreeMap<Long,Long> ref = new TreeMap(chr.getEncodeMap());
        //TreeMap<Long,Long> test = new TreeMap(read.getEncodeMap());
        SortedSet<Long> vals = new TreeSet(read.getEncodeMap().values()) ;
        int count = 0;
        for (Long val : vals){
            //System.out.println(key);
            System.out.println("Number of mapping iteration : " + count++);
//            System.out.println("the value is " + val + " Key is " + getKeyFromValue(read.getEncodeMap(),val));
            if (chr.getEncodeMap().containsKey(getKeyFromValue(read.getEncodeMap(),val))){
                System.out.println("Key " + getKeyFromValue(read.getEncodeMap(),val) +": Match at position " + chr.getEncodeMap().get(getKeyFromValue(read.getEncodeMap(),val)));
            }else{
                System.out.println(" Key " + getKeyFromValue(read.getEncodeMap(),val) + ": Not match");
            }   
            
        }
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
        return chr.getEncodeMap();
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

    
    
}
