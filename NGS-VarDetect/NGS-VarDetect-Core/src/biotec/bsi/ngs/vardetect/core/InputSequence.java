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
    
   private Vector<ShortgunSequence> seqs;
   private String name;
   
   public InputSequence(){
       seqs = new Vector<ShortgunSequence>();
   }
   
   public void addchrName(String chrName){
       this.name = chrName;
   }
   
   public void addRead(ShortgunSequence in){
       this.seqs.add(in);
   }
   
   public void addListOfRead(Vector<ShortgunSequence> inVector){
       this.seqs.addAll(inVector);
   }
   public String getChrName(){
       return this.name;
   }
   
   public Vector getInputSequence(){
       return this.seqs;
   }
   
   public int getNumberShortgunSequence(){
       return this.seqs.size();
   }
}
