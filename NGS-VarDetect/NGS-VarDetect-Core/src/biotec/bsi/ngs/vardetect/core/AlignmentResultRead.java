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
public class AlignmentResultRead {
    private ArrayList<ShortgunSequence> shrtRead;
    
    public AlignmentResultRead(){
        
        this.shrtRead = new ArrayList();
       
    }
    
    public void addResult(ShortgunSequence inRead){
        this.shrtRead.add(inRead);
    }
    
    public ArrayList<ShortgunSequence> getResult(){
    
        return this.shrtRead;
    }
        
    
}
