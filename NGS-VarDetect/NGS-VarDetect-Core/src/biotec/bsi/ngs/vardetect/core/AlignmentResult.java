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
    
    public void addResult(long idx, long chrNumber, String readName){
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
            
        
        //result.put(readName, code)
        
    }
}
