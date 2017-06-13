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
 */

public class TestBedCMD {
    
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
      
//       ReferenceSequence ref = SequenceUtil.readAndIndexReferenceSequence("/Users/soup/Desktop/hg19/hg19.fa");
        
//        String savefilename = "_SimulateLongRead";
//        System.out.println("Get reference sequence");
//        ReferenceSequence ref = SequenceUtil.getReferenceSequence(args[0],18); //runFile hg19.fa
//        //^^^^
//        //ChromosomeSequence c = ref.getChromosomeSequenceByName("chr21");
//        System.out.println("Simulate Data");
//        //InputSequence input =  SimulatorUtil_WholeGene.simulateWholeGene(ref, 5, 100, "20", "21");
//        //InputSequence input =  SimulatorUtil_WholeGene.simulateComplexWholeGeneRandom(ref,1, 100, 1);
//        InputSequence input =  SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomMixed(ref, 12, 30000, 2, 30000,100,savefilename);
        
//        Aligner aligner = AlignerFactory.getAligner();          // Will link to BinaryAligner
//          
//        AlignmentResultRead align = aligner.alignV3(ref, input,18,5);  // function align is located in binary aligner
//                
//        System.out.println("Do sortCountCutResult");
//        align.sortCountCutResultForMapV3(5);
//        System.out.println("Do write Report");
//        align.writeSortedCutResultMapToPathInFormatV3(ref.getPath(),savefilename, "txt");
//        System.out.println("Done");
//    
//        AlignmentResultRead readAlign = SequenceUtil.readAlignmentReportV2("/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/hg19"+savefilename+".txt",18);
//        System.out.println("Begin create color array");
//        Clustering.createColorArray(readAlign, 100, 18);
//        System.out.println("Done create color array");

          String fastaFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_050417/dm6_O2_4thread_th5_alignmentResult_longRead_Sample.fa";
          int startPoint = 0;
          int length = 15;
          SequenceUtil.truncateFastaFIles(fastaFile, startPoint, length);
          
//          SequenceUtil.miRNASeparator(fastaFile);
    }

}
