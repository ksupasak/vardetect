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
public class ReferenceSequence {

    String filename;
    
    public void setFilename(String filename) {

        this.filename = filename;
        
    }
    
    
    public String toString(){
        return "REF : "+this.filename;
    }
    
}
