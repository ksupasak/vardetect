/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;
import biotec.bsi.ngs.vardetect.core.*;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.FileNotFoundException;
import java.util.Enumeration;
import java.util.Vector;
/**
 *
 * @author soup
 */
public class NGSCMD {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException {
        // TODO code application logic here
        
       ReferenceSequence ref = SequenceUtil.readReferenceSequence(args[1]);
       
       Vector<ChromosomeSequence> chrs = ref.getChromosomes();
       
       Enumeration<ChromosomeSequence> e = chrs.elements();
       
       while(e.hasMoreElements()){
           
           ChromosomeSequence chr = e.nextElement();
           
           System.out.println(chr.getName());
            
           
//           EncodedSequence encode = SequenceUtil.encodeSerialChromosomeSequence(chr);
//           encode.writeToPath("/Users/soup/Desktop/hg19/"+chr.getName()+".map", "map");
           
           EncodedSequence encode = SequenceUtil.getEncodeSequence(chr);
           
           

       
       }
       
       
//        1111 = 15
          
//            System.out.println(""+((15*2)&15));
          



       
//         SequenceUtil.extractReferenceSequence(args[1], args[3]);
       
    }
    
}
