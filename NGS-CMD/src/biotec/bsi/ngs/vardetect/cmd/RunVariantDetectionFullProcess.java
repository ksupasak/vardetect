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
import biotec.bsi.ngs.vardetect.core.VariationResult;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 *
 * @author worawich
 */
public class RunVariantDetectionFullProcess {
    /**
     * It run cut repeat version of alignment (Clone from version 5 alignment function) But use alignment function that support long read alignment 
     * It use cut repeat protocol alignment function
     */
    
    public static void main(String[] args) throws IOException, InterruptedException {
        // TODO code application logic here
        /**
         * Alignment Process (Phase 1)
         */
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
        
        int overlap = Integer.valueOf(args[9]);
        byte percentMatch = Byte.valueOf(args[10]);
        int coverageThreshold = Integer.valueOf(args[11]); 
//        boolean annotationFlag = Boolean.valueOf(args[9]);              // Anotation flag true = do annotate ; false = not do
        
        System.out.println("Get reference sequence");
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer); //runFile hg19.fa
        
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
            InputSequence input = SequenceUtil.readSampleFileV4(inputPath,i,Math.min(numSample, i+propotion));
            
            Aligner aligner = AlignerFactory.getAligner();          // Will link to BinaryAligner

            AlignmentResultRead align = aligner.alignMultithreadLongReadRepeatCut(ref, input, numThread, numMer, threshold, repeatThreshold);  // function align is located in binary aligner
    
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
        
        /***********************************************************************/
        /**
         * Post  process (Phase 2)
         */
        String filename2 = folder.toString()+File.separator+filename;
//        String indexFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/drosophila/d.melanogaster/dm6_filter.index";            // use for traceback to natural chromosome name
//        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_050417/O3_S3_L001_R2_001.fa";
        String saveFileType = "txt";
        int numPart = count; // number of part of alignment result file
//        int readLength = 24;
        int merLength = numMer;
        int maxFullMatch = 100; // it's percentage value
//        int propotion = 1000000;
//        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/simLongRead/";
        
        
        File mainFile = new File(filename2);
        String path = mainFile.getParent();
        String saveFileName = mainFile.getName().split("part")[0]+"forLinuxSort";
        String sortedFileName = mainFile.getName().split("part")[0]+"Sort";
        String saveSampleFileName = mainFile.getName().split("part")[0]+"Sample";
//        File pathVar = new File(args[1]);
//        filename = pathVar.getName().split("\\.")[0];
//        path = pathVar.getParent()+"/";
//        saveFileName = pathVar.getPath().split("\\.")[0];
           
        for(int i=1;i<=numPart;i++){
            
            AlignmentResultRead readAlign = SequenceUtil.readAlignmentReportV2(filename2+i+".txt",merLength);
//            AlignmentResultRead readAlign = SequenceUtil.readBinaryAlignmentReportV3(filename+i+".bin",merLength);
            
            /**
             * Generate Bed graph  
             */
//            readAlign.countAlignMatch();
//            readAlign.writeMatchCountReport(filename+i);              // write bed graph file format for overall vistualize
//            readAlign.writeMatchCountReport(filename+i,indexFile);              // write bed graph file format for overall vistualize + index file
            /***************************************/
            
            System.out.println("Begin create color array");
            Clustering.createColorArrayV2(readAlign, merLength);        
            System.out.println("Done create color array");
//            readAlign.sortFilterResult("all",maxFullMatch);
            readAlign.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName, saveFileType,"all",maxFullMatch);
//            readAlign.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName, saveFileType, "all");
//*****for create set of sample that align********//            
//            readAlign.writeAlignSequenceReadFasta(path, saveSampleFileName, sampleFile, propotion, "all",maxFullMatch);
//*****************//
//            VariationResult varRes = SequenceUtil.analysisResultFromFile(path+saveFileNameForPostProcess+".txt",18,100);
            readAlign = null;
            System.gc();
        }
        
        /**
         * Implement execute sort command with linux command
         */
        String fullPathSaveUnSortFile = path+File.separator+saveFileName+"."+saveFileType;
        String fullPathSaveSortFile = path+File.separator+sortedFileName+"."+saveFileType;
        
        Process linuxControlProcess;
        try{
            String[] cmd = {"/bin/sh","-c","sort -t, -k10,10 -k9,9n "+fullPathSaveUnSortFile+" >> "+fullPathSaveSortFile};
            linuxControlProcess = Runtime.getRuntime().exec(cmd);
            linuxControlProcess.waitFor();
            String[] cmd2 = {"/bin/sh","-c","rm "+fullPathSaveUnSortFile};
            linuxControlProcess = Runtime.getRuntime().exec(cmd2);
            linuxControlProcess.waitFor();
        } catch (Exception e){
            e.printStackTrace();
        }
        
        /***********************************************************************/
        
        /**
         * Detect Variation (Phase 3)
         */
        String gffFile = "/Volumes/PromisePegasus/worawich/Referense/hg38/gff3/Homo_sapiens.GRCh38.87.chr.gff3";
//        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/hg38_sim_3/";
//        String saveFilename = "hg38_FullNewMethod_Sim_alignmentResult_VariantReport";
        
//        int readLength = 24;
//        int overlap = Integer.valueOf(args[9]);
//        byte percentMatch = Byte.valueOf(args[10]);
//        int coverageThreshold = Integer.valueOf(args[11]);       
//        String saveFilenameCov = filename + "_VariantCoverageReport_match" + percentMatch;
        
        VariationResult varRes = SequenceUtil.analysisResultFromFileV3(fullPathSaveSortFile,merLength,overlap,percentMatch);
        varRes.createVariantReport();
 
        varRes.analyzeCoverageFusionV2();
        varRes.sortCoverageFusion();
//        varRes.writeVariantCoverageVirtualizeWithAnnotationReportToFile(filename2, gffFile, coverageThreshold, 'F');
        varRes.writeVariantSortedCoverageVirtualizeReportToFile(filename2, coverageThreshold, 'F');
        varRes.writeVariantSortedCoverageReportToFile(filename2, coverageThreshold, 'F',true);
        
        varRes.analyzeCoverageIndelV2();
        varRes.sortCoverageIndel();
//        varRes.writeVariantCoverageVirtualizeWithAnnotationReportToFile(filename2, gffFile, coverageThreshold, 'I');
        varRes.writeVariantSortedCoverageVirtualizeReportToFile(filename2, coverageThreshold, 'I');
        varRes.writeVariantSortedCoverageReportToFile(filename2, coverageThreshold, 'I',true);

        System.gc();
        /***********************************************************************/
    }

    
}
