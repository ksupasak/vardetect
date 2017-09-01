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
import biotec.bsi.ngs.vardetect.core.util.VisualizeResult;
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
    
        String refPath = args[0];
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,18); //runFile hg19.fa
       
        InputSequence tempInSS = new InputSequence();
        int numMer = 18;
        String s = "GCTGGGATTACAGGCGTGAGCCACCGAGCCTGGCCAAACCATCACTTTTCATGAGCAGGGATGCACCCACTGGCACTCCTGCACCTCCCACCCTCCCCCT";
        
        for(int i=0;i<(s.length()-numMer)+1;i++){                                  // (Windowing with one stepping) for loop over String sequence which has limit round at (string length - mer length) + one [maximum possible mer sequence]
                int index = i;
                String sub = s.substring(i, i+numMer);                                 // cut String sequence into sub string sequence (mer length long) 
                //System.out.println("check sub length"+sub.length());
                long m = SequenceUtil.encodeMer(sub, numMer);
                
                System.out.println("index: "+i+" sequence: "+sub+" sequence code: "+m);
        }
        
        ShortgunSequence inSS = new ShortgunSequence(s);
        inSS.addReadName("err01");
        tempInSS.addRead(inSS);

//        tempInSS = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomMixed(ref, 200, 100, 10, 10000, 9);
        System.out.println("done");
                    
        Aligner aligner = AlignerFactory.getAligner();          // Will link to BinaryAligner

        AlignmentResultRead align = aligner.alignV3(ref, tempInSS,numMer,5);  // function align is located in binary aligner
        
        
    }
    
}