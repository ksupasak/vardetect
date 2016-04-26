/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

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
    
    
}
