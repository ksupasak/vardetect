/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

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
                
                ChromosomeSequence c = new ChromosomeSequence(chr,seq);
                
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

        ChromosomeSequence c = new ChromosomeSequence(chr,seq);
                
        ref.addChromosomeSequence(c);
        
    }
    
    
    
    
    } catch (IOException x) {
        System.err.format("IOException: %s%n", x);
    }    
           
           
     
       return ref;
   }
    





 public static void extractReferenceSequence(String filename, String topath){
     
       
       
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
                
                ChromosomeSequence c = new ChromosomeSequence(chr,seq);
                
                c.writeToPath(topath+chr,"FA");
                
                
                
            }
            seq = new StringBuffer();
            chr = line.substring(1,line.length());
            
        }else{
            
            seq.append(line.trim());
            
            
        }
        
    }
    
    if(seq.length()>0){
        
//        System.out.println("CHR : "+chr+" Size : "+seq.length());

        ChromosomeSequence c = new ChromosomeSequence(chr,seq);
       
        c.writeToPath(topath+chr,"FA");
        
    }
    
    
    
    
    } catch (IOException x) {
        System.err.format("IOException: %s%n", x);
    }    
           
           
     
   }


}