/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.alignment;

import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResult;
import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import java.io.IOException;
import java.util.Enumeration;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author soup
 */
public class BinaryAligner implements Aligner{

    ReferenceSequence ref;
    
    
    public void setReferenceSequence(ReferenceSequence ref){
        this.ref = ref;
        
    }
    

   
    public AlignmentResult align(ReferenceSequence ref, InputSequence input) {
        
        this.setReferenceSequence(ref);
        return align(input);
        
    }
    
    public AlignmentResult align(InputSequence input){
        
        
        AlignmentResult res = new AlignmentResult(input);
        
        int mer = 18;
        
        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){
            try {
                
                
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                EncodedSequence encoded = encodeSerialChromosomeSequenceV3(chr);
                
                while(seqs.hasMoreElements()){
                    ShortgunSequence seq = seqs.nextElement();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();
                    
                     System.out.print(chr.getName()+"\t");
                    
                    
                    for(int i=0;i<s.length()-mer;i++){
                        String sub = s.substring(i, i+mer);
                        long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){
                            m = m<<28;
                            long pos = encoded.align(m);
                            int idx = (int) (pos-i);
                            if(pos<0){
                              idx = 0;
                            }
//                            System.out.println(""+chr.getName()+" "+sub+" "+sub.length()+" : "+m+" pos : "+pos+" idx : "+idx);
                            
                            System.out.print("\t"+idx);
                            
                        }
                    
                    }
                     System.out.println();
                    
                    
                }
                System.gc();
                
                
                
            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
          
        
        
        return res;
    }
    
}
