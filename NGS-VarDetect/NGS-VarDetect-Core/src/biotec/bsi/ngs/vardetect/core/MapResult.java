/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

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
    
    
    public MapResult(){
    
        this.result = new HashMap(); 
    }
    
    public void addResult(Long inPos, Long numCount){
        
        chrnumber = inPos&255;
        alignPos = inPos>>8;
        
        if (result.containsKey(inPos)){
            numCount+=(result.get(inPos));
        }
        
        result.put(inPos,numCount);
    }
    
    public Map getResult(){
        return result;
    }
    
    
    
    
    
    
}
