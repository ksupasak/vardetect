/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Enumeration;
import java.util.List;
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
    
   
    public String getPath(){
//        String[] p = filename.split(File.pathSeparator);
//        String path = String.join(File.pathSeparator, p);
//        
        Path p = Paths.get(filename);
        Path folder = p.getParent();
        return folder.toString();
        
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
