/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.alignment.AlignerFactory;
import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ClusterGroup;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_WholeGene;
import biotec.bsi.ngs.vardetect.core.util.VisualizeResult;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class NGSCMD9 {
    public static void main(String[] args) throws IOException {
        
        System.out.println("Get reference sequence");
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(args[0],18); //runFile hg19.fa
        
        //ChromosomeSequence c = ref.getChromosomeSequenceByName("chr21");
        System.out.println("Read Data");
        //InputSequence input =  SimulatorUtil_WholeGene.simulateWholeGene(ref, 5, 100, "20", "21");
//        ShortgunSequence tempSS = SequenceUtil.readInputFile(args[1]);
       
        String sequence = "ATGGCACATGCAGCGCAAGTAGGTCTACAAGACGCTACTTCCCCTATCATAGAAGAGCTTATCACCTTTCATGAT" +
"CACGCCCTCATAATCATTTTCCTTATCTGCTTCCTAGTCCTGTATGCCCTTTTCCTAACACTCACAACAAAACTA" +
"ACTAATACTAACATCTCAGACGCTCAGGAAATAGAAACCGTCTGAACTATCCTGCCCGCCATCATCCTAGTCCTC" +
"ATCGCCCTCCCATCCCTACGCATCCTTTACATAACAGACGAGGTCAACGATCCCTCCCTTACCATCAAATCAATT" +
"GGCCACCAATGGTACTGAACCTACGAGTACACCGACTACGGCGGACTAATCTTCAACTCCTACATACTTCCCCCA" +
"TTATTCCTAGAACCAGGCGACCTGCGACTCCTTGACGTTGACAATCGAGTAGTACTCCCGATTGAAGCCCCCATT" +
"CGTATAATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAACAGATGCAATT" +
"CCCGGACGTCTAAACCAAACCACTTTCACCGCTACACGACCGGGGGTATACTACGGTCAATGCTCTGAAATCTGT" +
"GGAGCAAACCACAGTTTCATGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTT" +
"ACCCTATAG"; // got this Sequence from http://molevol.lysine.umiacs.umd.edu/resources/fileformats on the fasta unAlign Example file.
        ShortgunSequence inSS = new ShortgunSequence(sequence);
        inSS.addReadName("Un_Align_EX");
        InputSequence input = new InputSequence();
//        input.addRead(inSS);
        //input = SequenceUtil.readSampleFile(args[1]);
        String fixPath = args[1]; // input file path
        //String fixPath = "/Users/worawich/VMdev/3661/output.fa";
        input = SequenceUtil.readSampleFileV2(fixPath);
//        InputSequence input =  SimulatorUtil_WholeGene.simulateComplexWholeGeneRandom(ref,5, 100, 5);
        
        System.out.println("********** Do Alignment *********");
        Aligner aligner = AlignerFactory.getAligner();          // Will link to BinaryAligner
          
        AlignmentResultRead align = aligner.align(ref, input,18,5);  // function align is located in binary aligner
        
        
        AlignmentResultRead cutAlign = align.generateSortedCutResult(5);
        
//        Map<String,ArrayList<Map>> result = new HashMap();
//        result = align.getAlignmentResultV2();
        //ArrayList test = new ArrayList();
        //test = result.get("Read0");
        //System.out.print("/n");
        //System.out.print("Test represent Result: " + test.size());
        
//        VisualizeResult.visualizeAlignmentResultV2(align);
        
        //System.out.println("Size of Result: " + align.getAlignmentCount().size());
        
//        VisualizeResult.visualizeAlignmentCountMatchCutPlusColor(align,100);
        System.out.println("********** Do Write Result SortedResult *********");
        cutAlign.writeSortedResultToPath(ref.getPath(), "txt");
        System.out.println("********** Do Write Result UnSortedResult *********");
        cutAlign.writeUnSortedResultToPath(ref.getPath(), "txt");
        System.out.println("********** Do Write Result CutResult *********");
        cutAlign.writeSortedCutResultToPath(ref.getPath(), "txt", 5);
        
        cutAlign.writeSortedCutResultToPathInFormat(ref.getPath(), "txt", 5);
        
        AlignmentResultRead readAlign = SequenceUtil.readAlignmentReport("/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/hg19_Format_AlignSortedCutResult.txt");
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
        System.out.println("********** Do Calculate Euclidient Distance *********");
        cutAlign.calculateEuclidientdistance(); // Must have this order before clustering
//        VisualizeResult.visualizeDistanceTable(align);
        System.out.println("********** Do Write Result DistanceTable *********");
        cutAlign.writeDistanceTableToPath(ref.getPath(), "txt");
        
        System.out.println("********** Do Clustering Group *********");
        System.out.println("Number of Read in consider : "+ cutAlign.getResult().size());
        ArrayList<ClusterGroup> groupResult = Clustering.clusteringGroup(cutAlign, 100);
        cutAlign.addGroupReult(groupResult);
        System.out.println("********** Do Write Group Result *********");
        cutAlign.writeClusterGroupToPath(ref.getPath(), "txt");
        
        System.out.println(" check number of group : " + groupResult.size());
        VisualizeResult.visualizeClusterGoup(groupResult);
        
        System.out.println("********** Do Reconstruct Sequence *********");
        cutAlign.enableReconstruct();
        System.out.println("********** Do Write Pattern Report *********");
        cutAlign.writePatternReport(ref.getPath(), "txt");
        
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
