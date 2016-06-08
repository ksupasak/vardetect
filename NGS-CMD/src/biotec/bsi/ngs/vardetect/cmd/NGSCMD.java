/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;
import biotec.bsi.ngs.vardetect.alignment.AlignerFactory;
import biotec.bsi.ngs.vardetect.core.*;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
/**
 *
 * @author soup
 */
public class NGSCMD {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
//       ReferenceSequence ref = SequenceUtil.readAndIndexReferenceSequence("/Users/soup/Desktop/hg19/hg19.fa");
       
          ReferenceSequence ref = SequenceUtil.getReferenceSequence("/Users/soup/Desktop/hg19/hg19.fa");
          
          ChromosomeSequence c = ref.getChromosomeSequenceByName("chr21");

          InputSequence input =  SimulatorUtil.simulateIndel(c, 5, 100);
       
          
//          System.out.println(((1<<36)-1)<<28);
      // alignment
      
      
          Aligner aligner = AlignerFactory.getAligner();
          
          AlignmentResult align = aligner.align(ref, input);
      
      
//        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();
//
//        while(chrs.hasMoreElements()){
//            ChromosomeSequence chr = chrs.nextElement();
//            Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
//            EncodedSequence encoded = encodeSerialChromosomeSequenceV3(chr);
//            while(seqs.hasMoreElements()){
//                ShortgunSequence seq = seqs.nextElement();
//                System.out.println(""+chr.getName()+" ");
//            }
//            System.gc();
//            
//        }
//          
      //ReferenceSequence ref = SequenceUtil.readReferenceSequence(args[1]);
      
      
      //ChromosomeSequence chr = ref.getChromosomes().elementAt(0);
      //EncodedSequence encode = SequenceUtil.getEncodeSequence(chr);
      //System.out.println(chr.getName());
      //System.out.println(encode.getEncodeChrName());
      //String[] ss = encode.getEncodeChrName().split("r");
      //System.out.println(encode.getEncodeChrName().split("r"));
      
      //InputSequence is = SimulatorUtil.simulateIndel(chr, 5, 100);
      
      //Enumeration<ShortgunSequence> e = is.seqs.elements();
      
      //while(e.hasMoreElements()){
          
          //ShortgunSequence ss = e.nextElement();
          
          //EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(ss.seq);
          
          //SequenceUtil.mapGenome(encode, encodeSim);
          
          
          
          
      }
      
      
      
      
      
      
     
//      ReferenceSequence refB = SequenceUtil.readReferenceSequence(args[2]);
//       
//      Vector<ChromosomeSequence> chrs = refB.getChromosomes();
//       
//      Enumeration<ChromosomeSequence> e = chrs.ele&ments();
//      
//      //aonSystem.out.println(e.nextElement().getSequence().length());
//     
//      
//      while(e.hasMoreElements()){
//          
//          ChromosomeSequence chr = e.nextElement();
//          //ChromosomeSequence chr = e.nextElement();
//          
//          CharSequence test = SequenceUtil.concatenateChromosome(chr, chr, 100, 100);
//           
//          System.out.println(chr.getName());
//          System.out.println("Whole genome sim Data ="+test);
//            
//           
////           EncodedSequence encode = SequenceUtil.encodeSerialChromosomeSequence(chr);
////           encode.writeToPath("/Users/soup/Desktop/hg19/"+chr.getName()+".map", "map");
//           
//          EncodedSequence encode = SequenceUtil.getEncodeSequence(chr);
//          EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(test);
//          
//          System.out.println(encodeSim.getEncodeMap());
//          
//          //Hashtable genome = encode.getEncodeMap();
//          //Hashtable read = encodeSim.getEncodeMap();
//          //TreeMap read = encodeSim.getEncodeMap();
//          Map read = encodeSim.getEncodeMap();
//          
//          
//          SequenceUtil.mapGenome(encode, encodeSim);
//          
//          // Can match hashtable but still have some confusion (Is their any way to sorted the Key ?)
//       }
//       
//       
       
//        1111 = 15
          
//            System.out.println(""+((15*2)&15));
          



       
//         SequenceUtil.extractReferenceSequence(args[1], args[3]);
       
    //}
    
}
