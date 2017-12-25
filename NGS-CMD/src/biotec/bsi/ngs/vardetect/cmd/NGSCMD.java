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
    public static void main(String[] args) throws IOException, FileNotFoundException, InterruptedException {
        
         String refPath = args[0];
         CombineReferenceSequence ref = SequenceUtil.getCombineReferenceSequence(refPath,16); //runFile hg19.fa
         
         ref.setMinimumPeakPattern(10, 5);

// for large batch with multi thread         
//         ref.setNumberOfThread(8);
//         ref.setTotalRead(86000000);
//         ref.setSkipRead(0);
//         ref.setMaximumDuplicatePattern(3);

// for large batch with multi thread         
//         ref.setNumberOfThread(4);
//         ref.setTotalRead(100000);
//         ref.setSkipRead(0);
//         ref.setMaximumDuplicatePattern(3);
         
// for small batch with multi thread         
         ref.setNumberOfThread(4);
         ref.setTotalRead(10000);
         ref.setSkipRead(0);
         ref.setMaximumDuplicatePattern(1);
         ref.setRandomAccess(true);
         

// for debug         
         ref.setNumberOfThread(1);
         ref.setTotalRead(1);
         ref.setSkipRead(3997);
         ref.setMaximumDuplicatePattern(1);
         ref.setMinimumPeakPattern(10, 2);
         ref.setRandomAccess(true);
         
         
//         ref.searchMer(0);
//         ref.setOutputFile("/Users/soup/human/RB_cancer/277T_sorted_unmap.all.out");
//         ref.runProfileSV("/Users/soup/human/RB_cancer/277T_sorted_unmap.bam");
         
//         
//         
//         int found = 0;
//         for(int i=0;i<1000000;i++){
//         int v = (int)(Math.random()*Integer.MAX_VALUE);
//         int pos[] = ref.searchMer(v);
//         if(pos!=null){
//             found ++;
//            System.out.println(""+(i+1)+". Search for : "+Integer.toBinaryString(v)+" Pos "+ pos.length);
//         }else{
////            System.out.println(""+(i+1)+". Search for : "+Integer.toBinaryString(v)+" Not found");  
//         }
//         
//         }
//         System.out.println("Found for : "+found); 
         
//         int a= -100;
//         long b= ((long)a<<32)+3;//(a<<32);
//        
//         System.out.println("A = "+Integer.toBinaryString(a));
//         System.out.println("B = "+Long.toBinaryString(b));
         
         
         
         
    }
    /**
     * @param args the command line arguments
     */
    public static void main2(String[] args) throws IOException {
    
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