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
public class InputSequence {
    
   public Vector<ShortgunSequence> seqs;
   private String name;
   
   public InputSequence(){
       seqs = new Vector<ShortgunSequence>();
   }
   
   public void addName(String chrName){
       this.name = chrName;
   }
   public String getChrName(){
       return this.name;
   }
   
}
