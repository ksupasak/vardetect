/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 *
 * @author soup
 */
public class EncodedSequence {

    Hashtable<Long,Long> map;
    
    public void setMap(Hashtable<Long, Long> map) {
        this.map = map;
    }
    
    public Hashtable getEncodeMap(){
        return this.map;
    }
    public void readFromPath(String file_path, String fa) throws FileNotFoundException {
        
        map = new Hashtable<Long,Long>();
        
        
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
        
        
    }
    
    
    public void writeToPath(String path, String fa) throws FileNotFoundException {

       
       PrintStream ps = new PrintStream(path+"."+fa);
       
       Enumeration<Long> e = map.keys();
       
       while(e.hasMoreElements()){
           Long mer = e.nextElement();
           Long pos = map.get(mer);
           ps.println(mer+"\t"+pos);
         
       }
               
               

    }
    
}
