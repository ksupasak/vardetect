/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.Enumeration;
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
    
    public ChromosomeSequence getChromosomeSequenceByName(String name){
        return null;
    }
    
    
    public Vector<ChromosomeSequence> getChromosomes(){
        return chrs;
    }
    
   
    
    
    public String toString(){
        StringBuffer sb = new StringBuffer();
        
        sb.append("REF : "+this.filename+" # CHR : "+chrs.size()+"\n");
        
        Enumeration<ChromosomeSequence> e = chrs.elements();
        
        while(e.hasMoreElements()){
            ChromosomeSequence chr = e.nextElement();
            sb.append(chr.name+" : "+chr.seq.length()+"\n");
        }
        
        
        
        return sb.toString();
    }
    
    
    
}
