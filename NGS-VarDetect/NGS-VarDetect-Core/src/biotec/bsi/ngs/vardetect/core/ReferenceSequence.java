/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.Vector;

/**
 *
 * @author soup
 */
public class ReferenceSequence {

    String filename;
    Vector<ChromosomeSequence> chrs;
    
    public ReferenceSequence(){
        chrs = new Vector<ChromosomeSequence>();
    }
    
    public void setFilename(String filename) {

        this.filename = filename;
        
    }
    
    public void addChromosomeSequence(ChromosomeSequence chr){
        chrs.add(chr);
    }
    
    
    public String toString(){
        return "REF : "+this.filename+" CHR : "+chrs.size();
    }
    
    
    
}
