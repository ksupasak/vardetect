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
 * @author soup
 */
public class AlignmentResult {
    
    Map<Long,long[]> merPosMap;
    Map<String,ArrayList<Map>> resultMap;
    ArrayList<Map> arrayMap;
            
    Map<String,Long> result;
    
    Map<String,ArrayList> newResult;
    private long index;
    private long chrNum;
    private long code = 0;
    private long count = 0;

    ArrayList tempResult;    
    ArrayList<Long> arrayResult;
    ArrayList chrNumber;
    ArrayList alignPosition;
    ArrayList numMatch;
    ArrayList readName;

    InputSequence input;
    
    public  AlignmentResult(InputSequence input){
        this.merPosMap = new HashMap();
        this.resultMap = new HashMap();
        this.arrayMap = new ArrayList();
        
        this.input = input;  
        this.result = new HashMap();
        
        
        this.newResult = new HashMap();
        this.tempResult = new ArrayList();         
        this.arrayResult = new ArrayList();
        this.alignPosition = new ArrayList();
        this.chrNumber = new ArrayList();
        this.numMatch = new ArrayList();
        this.readName = new ArrayList();    
    }
    
    public Map addResult(long mer, long chrNumber, long[] pos){
        
        int len;
        long[] code = pos; 
        
        if(pos!=null&&pos.length>0){
            len = pos.length;
            if(pos[0] > 0){
                for(int i=0;i<len;i++){
                    code[i] = (chrNumber<<28)+pos[i];
                }
                System.out.println(" this is code check: " + code[0]);
                System.out.println(" this is mer check: " + mer);
                System.out.println(merPosMap == null);
                this.merPosMap.put(mer, code);
            }
        }
        
        return this.merPosMap;
        
        /*
        this.newResult.put(readName, this.arrayResult);
        this.index = idx;
        long dummy_subcode = (idx<<8)+chrNumber;
        
        if (dummy_subcode != this.code){
           
            if (count != 0){
                
                this.arrayResult = this.newResult.get(readName);
                this.code = (dummy_subcode<<16)+count;
                this.arrayResult.add(this.code);
                
                this.newResult.put(readName,this.arrayResult);
            }
            ///// Problem: how to collect the last result because we check variable changing but there is no change for the last result
            
            this.code = dummy_subcode;
            count = 1;
            
            //tempResult.
            // Create Array to collect pre result in case that one read has more than one match
        }
        else if(dummy_subcode == this.code){
            count++;
        }
            
        */
        //result.put(readName, code)
        
    }
    
    public void createMap(String readName, Map merMap){
        if (merMap != null){
            this.merPosMap = new HashMap(); // reset map do it every time 

            if (resultMap.containsKey(readName)){
                ArrayList<Map> dummy_merMap = new ArrayList();
                dummy_merMap = resultMap.get(readName);
                dummy_merMap.add(merMap);
                resultMap.put(readName, dummy_merMap);
            }
            else{
                ArrayList<Map> dummy_merMap = new ArrayList();
                dummy_merMap.add(merMap);                 
                resultMap.put(readName,dummy_merMap);
            }  
        }
    }
    
    public Map getAlignmentResult(){
        return this.resultMap;
    }
    
    
}
