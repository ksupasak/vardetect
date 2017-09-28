/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.alignment.AlignerFactory;
import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_WholeGene;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 *
 * @author worawich
 */
public class RunEvaluateFullProcessVariantDetection {
    /**
     * This function run full process of variant detection algorithm
     * 1. import input (import file or simulate data) simulate data is use for evaluate algorithm
     * 2. alignment
     * 3. post processing
     * 4. variant detection
     * 5. evaluate  result
     * 
     * It run cut repeat version of alignment (Clone from version 5 alignment function) But use alignment function that support long read alignment 
     * It use cut repeat protocol alignment function
     * 
     * 
     */
    
    public static void main(String[] args) throws IOException, InterruptedException {
        
        long startAlignTime = System.currentTimeMillis();
        String refPath = args[0];                                       // First argument; indicate reference  file (include it path if it not in the current directory)
        String inputPath = args[1];                                     // Second argument; indicate input file (include it path if it not in the current directory)
        String filename = args[2];                                      // Third argument; indicate save file name
        int propotion = Integer.valueOf(args[3]);                       // Forth argument; indicate the number of read per time
        int numMer = Integer.valueOf(args[4]);
        int threshold = Integer.valueOf(args[5]);                       // Fifth argument; indicate count number of mer match threshold
        int numThread = Integer.valueOf(args[6]);                       // Sixth argument; indicate number of thread
        int repeatThreshold = Integer.valueOf(args[7]);                 // Seventh argument; indicate number of repeat that we can accept. If the number of repeat is more than this threshold it will be cut out
                                                                   // If repeatThreshold has set to zero the protocol will change to not cut repeat protocol (consider all repeat)
        String filetype = args[8];
        
        int maxMultiplier = 5;              // use to be the exponent of indelSizeBase
        int indelSizeBase = 10;     
        char variantType = 'L';             // 4 variant type F=fusion L=large indel I=insert D=delete
        char insertSNPFlag = 'T';               // T = true(insert) and F = false(not insert)
        int numRead = 20;
        int readLen = 100;
        int readCoverage = 30;
        int diffL = 10000;
        int indelSize = 10;
        
        for(int indelSizeMultiplier=1; indelSizeMultiplier<maxMultiplier;indelSizeMultiplier++){        // indelSizeMultiplier start at 1
            long startTime = System.currentTimeMillis();
            
            int minIndelSize = (int)Math.pow(indelSizeBase, indelSizeMultiplier);                       // minIndelSize is indelSizeBase power by indelSizeMultiplier
            int maxIndelSize = (int)Math.pow(indelSizeBase, indelSizeMultiplier+1);                     // maxIndelSize is indelSizeBase power by indelSizeMultiplier+1
            
            
            String sampleFilename = "hg38_sim_"+variantType+indelSizeMultiplier; 
            ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer); //runFile reference
            /**
             * Simulate data
             */
            InputSequence input = new InputSequence();
    //        tempInSS = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomSingleType(ref, numRead, readLen, readCoverage, diffL, indelSize,filename,variantType);
            input = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomSingleTypeFixRange(ref, numRead, readLen, readCoverage, minIndelSize, maxIndelSize, diffL, indelSize,sampleFilename,variantType,insertSNPFlag);
            System.out.println("done");
            /*********************/
        }
//            /**
//             * Alignment
//             */
//
//            String savefilename = filename+"_indelSizeMultiplier"+indelSizeMultiplier;
//
//            Aligner aligner = AlignerFactory.getAligner();          // Will link to BinaryAligner
//
//            AlignmentResultRead align = aligner.alignMultithreadLongReadRepeatCut(ref, input, numThread, numMer, threshold, repeatThreshold);  // function align is located in binary aligner
//
//            long stopAlignTime = System.currentTimeMillis();
//            double totalAlignTime = ((stopAlignTime - startAlignTime)/1000)/60;
//            System.out.println(String.format("Alignment Time use : %.4f min",totalAlignTime));
//            System.out.println("Do write Report");            
////            align.writeSortedCutResultMapToPathInFormatLongRead(folder.toString()+File.separator,savefilename, filetype);
//            System.out.println("Done indelSizeMultiplier number " + indelSizeMultiplier);
//
//            long endTime = System.currentTimeMillis();
//            double totalTime = ((endTime - startTime)/1000)/60;
//            System.out.println(String.format("Time use : %.4f min",totalTime));
//            System.out.println();
//            /*********************/
//            
//            /**
//             * PostProcess
//             */
//            
//            String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/Thalasemia_Data/3662/3662_unmapped_alignResult_part";
//            String indexFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/drosophila/d.melanogaster/dm6_filter.index";            // use for traceback to natural chromosome name
//            String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_050417/O3_S3_L001_R2_001.fa";
//            String saveFileType = "txt";
//            int numPart = 2;
//    //        int readLength = 24;
//            int maxFullMatch = 100; // it's percentage value
//            
//            File mainFile = new File(filename);
//            String path = mainFile.getParent();
//            String saveFileName = mainFile.getName().split("part")[0]+"forLinuxSort";
//            
//
//
//                /**
//                 * Generate Bed graph  
//                 */
//    //            readAlign.countAlignMatch();
//    //            readAlign.writeMatchCountReport(filename+i);              // write bed graph file format for overall vistualize
//    //            readAlign.writeMatchCountReport(filename+i,indexFile);              // write bed graph file format for overall vistualize + index file
//                /***************************************/
//
//                System.out.println("Begin create color array");
//                Clustering.createColorArrayV2(align, numMer);        
//                System.out.println("Done create color array");
//    //            readAlign.sortFilterResult("all",maxFullMatch);
//                align.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName, saveFileType,"all",maxFullMatch);
//    //            readAlign.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName, saveFileType, "all");
//    //*****for create set of sample that align********//            
//    //            readAlign.writeAlignSequenceReadFasta(path, saveSampleFileName, sampleFile, propotion, "all",maxFullMatch);
//    //*****************//
//    //            VariationResult varRes = SequenceUtil.analysisResultFromFile(path+saveFileNameForPostProcess+".txt",18,100);
//                align = null;
//                System.gc();
//        }
//        }
//    
    }
}
