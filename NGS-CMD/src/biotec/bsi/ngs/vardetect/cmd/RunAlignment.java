/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.alignment.AlignerFactory;
import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResult;
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.ClusterGroup;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
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
 * 
 * Clone of NGSCMD 4 (Test bed for new align implement)
 */

public class RunAlignment {
    
     public static void main(String[] args) throws IOException {
        // TODO code application logic here
        String refPath = args[0];                                       // First argument; indicate reference  file (include it path if it not in the current directory)
        String inputPath = args[1];                                     // Second argument; indicate input file (include it path if it not in the current directory)
        int propotion = Integer.valueOf(args[2]);                       // Third argument; indicate the number of read per time
        
//       ReferenceSequence ref = SequenceUtil.readAndIndexReferenceSequence("/Users/soup/Desktop/hg19/hg19.fa");
        Map<String,ArrayList<Map>> aon = new HashMap();
        System.out.println(aon==null);
        System.out.println(aon!=null);
        
        System.out.println("Get reference sequence");
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,18); //runFile hg19.fa
        
        //ChromosomeSequence c = ref.getChromosomeSequenceByName("chr21");
        System.out.println("Simulate Data");
        //InputSequence input =  SimulatorUtil_WholeGene.simulateWholeGene(ref, 5, 100, "20", "21");
        //InputSequence input =  SimulatorUtil_WholeGene.simulateComplexWholeGeneRandom(ref,1, 100, 5);
        
        //InputSequence input = new InputSequence();
//        input.addRead(inSS);
        //input = SequenceUtil.readSampleFile(args[1]);
        
        String fixPath = "/Users/worawich/VMdev/3661/output.fa";
        int numSample = SequenceUtil.getNumberSample(inputPath);
        
        int count = 0;
        System.out.println("Total Sample: " + numSample);
        System.out.println("Propotion " + propotion + " read per part");
        for (int i = 0 ; i < numSample ; i += propotion){                       // loop over the input sample ( number of loop is up to the number of read per time )
            count++;
            String savefilename = "_Format_AlignSortedCutResultMap_part"+count;
            InputSequence input = SequenceUtil.readSampleFileV2(inputPath,i,Math.min(numSample, i+propotion));
            //input = SequenceUtil.readSampleFileV2(fixPath);


            Aligner aligner = AlignerFactory.getAligner();          // Will link to BinaryAligner

            AlignmentResultRead align = aligner.alignV2(ref, input,18,5);  // function align is located in binary aligner


    //        Map<String,ArrayList<Map>> result = new HashMap();
    //        result = align.getAlignmentResultV2();
            //ArrayList test = new ArrayList();
            //test = result.get("Read0");
            //System.out.print("/n");
            //System.out.print("Test represent Result: " + test.size());

    //        VisualizeResult.visualizeAlignmentResultV2(align);

            //System.out.println("Size of Result: " + align.getAlignmentCount().size());

    //        VisualizeResult.visualizeAlignmentCountMatchCutPlusColor(align,100);
            System.out.println("Do sortCountCutResult");
            align.sortCountCutResultForMapV2(5);
            System.out.println("Do write Report");
            align.writeSortedCutResultMapToPathInFormat(ref.getPath(),savefilename, "txt");
            System.out.println("Done part " + count);
        }
//        align.writeSortedResultToPath(ref.getPath(), "txt");
//        align.writeUnSortedResultToPath(ref.getPath(), "txt");
//        align.writeSortedCutResultToPath(ref.getPath(), "txt", 5);
        
//        cutAlign.writeSortedResultToPath(ref.getPath(), "txt");
//        cutAlign.writeUnSortedResultToPath(ref.getPath(), "txt");
//        cutAlign.writeSortedCutResultToPath(ref.getPath(), "txt", 5);
        
        /* Old Grouping algorithm
        align.createAllClusterCode();
        align.createAllClusterCodeSorted();
        
        for(int i =0;i<align.getAllClusterCode().length;i++){
            System.out.println("Check cluster code: " + align.getAllClusterCode()[i]);
        }
        for(int i =0;i<align.getAllClusterCode().length;i++){
            System.out.println("Check cluster code sorted: " + align.getAllClusterCodeSorted()[i]);
        }
        
        align.createGroupingResult();
        System.out.println(" ****** Check cluster result: " + align.getclusterResult().size());
        */
        
//        align.calculateEuclidientdistance(); // Must have this order before clustering
//        VisualizeResult.visualizeDistanceTable(align);
//        align.writeDistanceTableToPath(ref.getPath(), "txt");
//        
//        
//        ArrayList<ClusterGroup> groupResult = Clustering.clusteringGroup(align, 100);
//        align.addGroupReult(groupResult);
//        align.writeClusterGroupToPath(ref.getPath(), "txt");
//        
//        System.out.println(" check number of group : " + groupResult.size());
//        VisualizeResult.visualizeClusterGoup(groupResult);
//        
//        align.enableReconstruct();
//        align.writePatternReport(ref.getPath(), "txt");
        
        // Create save result path by just plugin AlignmentResult
        
//        System.out.println(" Size new resutl check " + align.getResult().size());
//        System.out.println(" Data check "+align.getResult().get(0).getReadName());
//        System.out.println(" Data check "+align.getResult().get(0).getSequence());
//        
//        System.out.println(" Data check align result map empty or not : "+align.getResult().get(0).getAlignmentCount().isEmpty());
        //VisualizeResult.visualizeAlignmentResult(align);
        
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
