/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class MerRead {
    
    private long merCode;
    private int index;
    private ArrayList<Long> chrPos;
    private ArrayList<Long> chrAlgn;
    
    public MerRead(){
        this.chrPos = new ArrayList();
        this.chrAlgn = new ArrayList();
    
    }
    
    // pos must be a compose number of chr:position and array of long (long[])
    public void addMatchResult(long mer,long[] pos, int idx){
        
        this.merCode = mer;
        this.index = idx;
        
        if(pos!=null&&pos.length>0){
            int len = pos.length;
            if(pos[0] > 0){
                for(int i=0;i<len;i++){
                    this.chrPos.add(pos[i]);
                }
            }
        }                
    }
    
}
