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
public class CombinedPos {
    

      
       String chr;
       int pos;
       
       CombinedPos(String chr, int pos){
           this.chr = chr;
           this.pos = pos;
       }
       
       public String getChr(){return this.chr;}
       public int getPos(){return this.pos;}
       
      
   }