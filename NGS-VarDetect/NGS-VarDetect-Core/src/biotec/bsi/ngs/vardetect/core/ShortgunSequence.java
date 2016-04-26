/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

/**
 *
 * @author soup
 */
public class ShortgunSequence {
 
    private String seq;
    private String readName;
    
    public ShortgunSequence(String seq){
        this.seq = seq;
    }
    
    public void addReadName(String readName){
        this.readName = readName;
    }
    
    public String getReadName(){
        return readName;
    }
    public String getSequence(){
        return seq;
    }
    
}
