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
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
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

public class RunAlignmentMultiThreadLongRead {
    /**
     * It run cut repeat version of alignment (Clone from version 5 alignment function) But use alignment function that support alignment long read 
     */
    
    public static void main(String[] args) throws IOException, InterruptedException {
        // TODO code application logic here
        long startAlignTime = System.currentTimeMillis();
        String refPath = args[0];                                       // First argument; indicate reference  file (include it path if it not in the current directory)
        String inputPath = args[1];                                     // Second argument; indicate input file (include it path if it not in the current directory)
        String filename = args[2];                                      // Third argument; indicate save file name
        int propotion = Integer.valueOf(args[3]);                       // Forth argument; indicate the number of read per time
        int numMer = Integer.valueOf(args[4]);
        int threshold = Integer.valueOf(args[5]);                       // Fifth argument; indicate count number threshold
        int numThread = Integer.valueOf(args[6]);                       // Sixth argument; indicate number of thread
        String filetype = args[7]; 
        
//       ReferenceSequence ref = SequenceUtil.readAndIndexReferenceSequence("/Users/soup/Desktop/hg19/hg19.fa");

        System.out.println("Get reference sequence");
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer); //runFile hg19.fa
        
        //ChromosomeSequence c = ref.getChromosomeSequenceByName("chr21");
        //System.out.println("Simulate Data");
        //InputSequence input =  SimulatorUtil_WholeGene.simulateWholeGene(ref, 5, 100, "20", "21");
        //InputSequence input =  SimulatorUtil_WholeGene.simulateComplexWholeGeneRandom(ref,1, 100, 5);
        
        //InputSequence input = new InputSequence();
//        input.addRead(inSS);
        //input = SequenceUtil.readSampleFile(args[1]);
        
        //String fixPath = "/Users/worawich/VMdev/3661/output.fa";
        Path inPath = Paths.get(inputPath);
        Path folder = inPath.getParent();
        int numSample = SequenceUtil.getNumberSample(inputPath);
        
        int count = 0;
        System.out.println("Total Sample: " + numSample);
        System.out.println("Propotion " + propotion + " read per part");
        for (int i = 0 ; i < numSample ; i += propotion){                       // loop over the input sample ( number of loop is up to the number of read per time )
            long startTime = System.currentTimeMillis();
            count++;
            String savefilename = filename+count;
            InputSequence input = SequenceUtil.readSampleFileV3(inputPath,i,Math.min(numSample, i+propotion));
            //input = SequenceUtil.readSampleFileV2(fixPath);


            Aligner aligner = AlignerFactory.getAligner();          // Will link to BinaryAligner

            AlignmentResultRead align = aligner.alignMultithreadLongRead(ref, input, numThread, numMer, threshold);  // function align is located in binary aligner


    
//            System.out.println("Do sortCountCutResult");
//            align.sortCountCutResultForMapV3(threshold);
            long stopAlignTime = System.currentTimeMillis();
            double totalAlignTime = ((stopAlignTime - startAlignTime)/1000)/60;
            System.out.println(String.format("Alignment Time use : %.4f min",totalAlignTime));
            System.out.println("Do write Report");
            align.writeSortedCutResultMapToPathInFormatLongRead(folder.toString()+File.separator,savefilename, filetype);
            System.out.println("Done part " + count);
            
            long endTime = System.currentTimeMillis();
            double totalTime = ((endTime - startTime)/1000)/60;
            System.out.println(String.format("Time use : %.4f min",totalTime));
            System.out.println();
        }
    
    }

}
