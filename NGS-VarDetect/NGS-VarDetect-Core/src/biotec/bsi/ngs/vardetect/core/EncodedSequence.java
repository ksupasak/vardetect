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
    Map<Long,Long> map;
    
    public void setMap(Map<Long, Long> map) {
        this.map = map;
    }
    
    public Map getEncodeMap(){
        return this.map;
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
        
        }else
        if(fa.compareTo("bmap")==0){
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
           
       }}
       else
       if(fa.compareTo("bmap")==0){
           
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
               
               

}
