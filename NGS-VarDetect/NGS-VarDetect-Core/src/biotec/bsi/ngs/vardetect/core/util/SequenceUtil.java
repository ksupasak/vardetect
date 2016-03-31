/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ReferenceSequence;

/**
 *
 * @author soup
 */
public class SequenceUtil {
   
   public static ReferenceSequence  readReferenceSequence(String filename){
       ReferenceSequence ref = new ReferenceSequence();
       ref.setFilename(filename);
       return ref;
   }
    
}
