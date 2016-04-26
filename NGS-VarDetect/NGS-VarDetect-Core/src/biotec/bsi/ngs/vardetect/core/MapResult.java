/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
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

    Map<Long,Long> keypos;
    Map<Long,Long> result;
    ArrayList arrayResult;
    ArrayList chrNumber;
    ArrayList alignPosition;
    ArrayList numMatch;
    ArrayList readName;
    
    
    public MapResult(){
    
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
        
        for(Map.Entry<Long,Long> entry : result.entrySet()){
            
            alignPosition.add(entry.getKey());
            numMatch.add(entry.getValue());
            chrNumber.add((entry.getKey())&255);
            readName.add(shotgunName);
        
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
    
    public void writeToPath(String path, String fa) throws FileNotFoundException, IOException {
        
        


        //Enumeration<Long> e = map.keys();
        if(fa.compareTo("map")==0){
            PrintStream ps = new PrintStream(path+"/Map_Result."+fa);
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
