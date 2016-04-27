/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author worawich
 */



public class MapResult {
    long alignPos;
    long key;
    long chrnumber;
    long count;
    String mapResultPath;

    Map<Long,Long> keypos;
    Map<Long,Long> result;
    Map<String,ArrayList> newResult;
    
    ArrayList inResult;    
    ArrayList arrayResult;
    ArrayList chrNumber;
    ArrayList alignPosition;
    ArrayList numMatch;
    ArrayList readName;
    
    
    public MapResult(){
        
        this.newResult = new HashMap();
        this.inResult = new ArrayList();
        
        this.result = new HashMap();
        this.arrayResult = new ArrayList();
        this.alignPosition = new ArrayList();
        this.chrNumber = new ArrayList();
        this.numMatch = new ArrayList();
        this.readName = new ArrayList();
    }
    
    public void addResultMap(Long inPos, Long numCount){
        
        chrnumber = inPos&255;
        alignPos = inPos>>8;
        
        if (result.containsKey(inPos)){
            numCount+=(result.get(inPos));
        }
        
        result.put(inPos,numCount);
    }
    
    public void addResultArray(Map inMap){
        
        arrayResult.add(inMap);
     
    }
    
    public void addResult(Map resultMap,String shotgunName){
        
        this.result = resultMap;
        String dummy;
        for(Map.Entry<Long,Long> entry : result.entrySet()){
            
            alignPosition.add(entry.getKey());
            numMatch.add(entry.getValue());
            chrNumber.add((entry.getKey())&255);
            readName.add(shotgunName);
            dummy = entry.getKey()+"\t"+((entry.getKey())&255)+"\t"+entry.getValue();     
                     
            
        }
    }
    public ArrayList getAlignPosition(){
        return alignPosition;
    }
    public ArrayList getNumMatch(){
        return numMatch;
    }
     public ArrayList getchrNumber(){
        return chrNumber;
    }
     public ArrayList getReadName(){
        return readName;
    }
    
    public Map getResultMap(){
        return result;
    }
    
    public ArrayList getResultArray(){
        return arrayResult;
    }
    
    public void getResult(){
        
        for (long i = 0;i<=arrayResult.size();i++){
             
        } 
    }
    public void readFromPath(String file_path,String fa){
        
        this.alignPosition = new ArrayList();
        this.chrNumber = new ArrayList();
        this.numMatch = new ArrayList();
        this.readName = new ArrayList();
        
        if(fa.compareTo("map")==0){
        
            Charset charset = Charset.forName("US-ASCII");

            Path path = Paths.get(file_path);


            StringBuffer seq = new StringBuffer();

            try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
                String line = null;
                int count = 0 ;



                while ((line = reader.readLine()) != null) {
                    String[] st = line.split("\t");
                    
                    String readN = st[0];
                    long alignAt = Long.valueOf(st[1]);
                    long chrN = Long.valueOf(st[2]);
                    long numM = Long.valueOf(st[3]);

                    alignPosition.add(alignAt);
                    chrNumber.add(chrN);
                    numMatch.add(numM);
                    readName.add(readN);
                    
                   

                    if(count%1000000==0)System.out.println("Read Mer "+count);
                    count ++;
                }   
            //System.out.println("Total mer : "+map.size());

            }
            catch(Exception e){
                    
            }
        }
    }
    
    public void writeToPath(String path, String fa) throws FileNotFoundException, IOException {
        
        

        this.mapResultPath = path+"/Map_Result."+fa;
        //Enumeration<Long> e = map.keys();
        if(fa.compareTo("map")==0){
            PrintStream ps = new PrintStream(mapResultPath);
            for (int i = 0;i<alignPosition.size();i++){
                ps.println(readName.get(i)+"\t"+alignPosition.get(i)+"\t"+chrNumber.get(i)+"\t"+numMatch.get(i));
                
            }
        }
        else if(fa.compareTo("bmap")==0){
        
            /*DataOutputStream os = new DataOutputStream(new FileOutputStream(path+"."+fa));
            System.out.println("Total line in file : " + alignPosition.size());
            for (int i = 0;i<alignPosition.size();i++){
              
              os.writeBytes(fa);
              System.out.println("Read name : " + mapRes.getReadName().get(i) + " Align at : " + mapRes.getAlignPosition().get(i) + " : On chromosome : "+ mapRes.getchrNumber().get(i)+ " : Match : " + mapRes.getNumMatch().get(i));

            }

            os.writeInt(map.keySet().size());
        
            for (Map.Entry<Long,Long> entry : map.entrySet()){
                Long mer = entry.getKey();
                Long pos = map.get(mer);
                os.writeLong(mer);
                os.writeLong(pos);
            }
        os.close();  */  

        }    
    } 
    
}
