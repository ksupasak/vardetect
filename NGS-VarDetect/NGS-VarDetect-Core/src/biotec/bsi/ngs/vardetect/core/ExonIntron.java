/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.Vector;

/**
 *
 * @author worawich
 */
public class ExonIntron {
    
    String chrName;
    String geneName;
    long startPos;
    long stopPos;
    int direction;
    
    //Vector<ChromosomeSequence> chrs ;
    public ExonIntron(String chrName,String geneName,long startPos,long stopPos,int direction){
        
        this.chrName = chrName;
        this.geneName = geneName;
        this.startPos = startPos;
        this.stopPos = stopPos;
        this.direction = direction;
    }
    
    public String getChrName(){
        return chrName;
    }
    public String getGeneName(){
        return geneName;
    }
    public long getStartPos(){
        return startPos;
    }
     public long getStopPos(){
        return stopPos;
    }
    public int getdirection(){
        return direction;
    }
    
}
