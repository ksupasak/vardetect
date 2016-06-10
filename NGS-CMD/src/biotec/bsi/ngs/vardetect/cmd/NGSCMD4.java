/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.alignment.AlignerFactory;
import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResult;
import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_WholeGene;
import biotec.bsi.ngs.vardetect.core.util.VisualizeResult;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */

public class NGSCMD4 {
    
     public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
//       ReferenceSequence ref = SequenceUtil.readAndIndexReferenceSequence("/Users/soup/Desktop/hg19/hg19.fa");
        Map<String,ArrayList<Map>> aon = new HashMap();
        System.out.println(aon==null);
        System.out.println(aon!=null);
        
        System.out.println("Get reference sequence");
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(args[0]);
       
        //ChromosomeSequence c = ref.getChromosomeSequenceByName("chr21");
        System.out.println("Simulate Data");
        InputSequence input =  SimulatorUtil_WholeGene.simulateWholeGene(ref, 5, 100, 20, 21);
        
        
        
        Aligner aligner = AlignerFactory.getAligner();
          
        AlignmentResult align = aligner.align(ref, input);
        
        
        Map<String,ArrayList<Map>> result = new HashMap();
        result = align.getAlignmentResult();
        ArrayList test = new ArrayList();
        test = result.get("Read0");
        System.out.print("/n");
        System.out.print("Test represent Result: " + test.size());
        
        VisualizeResult.visualizeAlignmentResult(align);
        
        // Pass!! Next create represent data part //
        
       /* ChromosomeSequence aon = ref.getChromosomeSequenceByName("chr21");
      // alignment
        EncodedSequence test = SequenceUtil.getEncodeSequenceV2(aon);
        System.out.println(test.getMers());*/
      
        /*Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){
            ChromosomeSequence chr = chrs.nextElement();
            Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
            EncodedSequence encoded = encodeSerialChromosomeSequenceV3(chr);
            while(seqs.hasMoreElements()){
                ShortgunSequence seq = seqs.nextElement();
                System.out.println(""+chr.getName()+" ");
            }
            System.gc();
            
        }*/
    
    }

}
