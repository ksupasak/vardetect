/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Enumeration;

/**
 *
 * @author Worawich
 */
public class NGSCMD3 {
    
    public static void main(String[] args){
        
      //ReferenceSequence ref = SequenceUtil.readReferenceSequence(args[1]);
      ReferenceSequence ref = SequenceUtil.readReferenceSequence(args[1]);
      
      
      ChromosomeSequence chr = ref.getChromosomes().elementAt(0);
      EncodedSequence encode = SequenceUtil.encodeSerialChromosomeSequence(chr);
      
      InputSequence is = SimulatorUtil.simulateIndel(chr, 5, 100);
      
      Enumeration<ShortgunSequence> e = is.seqs.elements();
      
      while(e.hasMoreElements()){
          
          ShortgunSequence ss = e.nextElement();
          
          EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(ss.seq);
          
          SequenceUtil.mapGenome(encode, encodeSim);
          
          
      }
        
    }
    
}
