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
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_WholeGene;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 *
 * @author worawich
 */
public class TestPostProcess {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
//        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/hg38_sim_L_SNP_F_5/hg38_sim_L_SNP_F_5_alignResult_part";
        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/SLE/SLE-9/SLE-9_hg38_part";
        String indexFile = "/Volumes/PromisePegasus/worawich/Referense/TB_reference_12Mer/H37Rv_NC_000962_reIndex.index";            // use for traceback to natural chromosome name
        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_050417/O3_S3_L001_R2_001.fa";
        String saveFileType = "txt";
        int numPart = 4;
//        int readLength = 24;
        int merLength = 18;
        int maxFullMatch = 100; // it's percentage value
        int propotion = 1000000;
//        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/simLongRead/";
        
        
        File mainFile = new File(filename);
        String path = mainFile.getParent();
        String saveFileName = mainFile.getName().split("part")[0]+"forLinuxSort";
        String sortedFileName = mainFile.getName().split("part")[0]+"Sort";
        String saveSampleFileName = mainFile.getName().split("part")[0]+"Sample";
//        File pathVar = new File(args[1]);
//        filename = pathVar.getName().split("\\.")[0];
//        path = pathVar.getParent()+"/";
//        saveFileName = pathVar.getPath().split("\\.")[0];
           
        for(int i=1;i<=numPart;i++){
            
            AlignmentResultRead readAlign = SequenceUtil.readAlignmentReportV2(filename+i+".txt",merLength);
//            AlignmentResultRead readAlign = SequenceUtil.readBinaryAlignmentReportV3(filename+i+".bin",merLength);
            
            /**
             * Generate Bed graph  
             */
//            readAlign.countAlignMatch();
//            readAlign.writeMatchCountReport(filename+i);              // write bed graph file format for overall vistualize
//            readAlign.writeMatchCountReport(filename+i,indexFile);              // write bed graph file format for overall vistualize + index file (indexfile should be filename.index not filename.fa.index)
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
        /*************************************************************************************/
        
        /**
         * Detect Variation
         */
        String gffFile = "/Volumes/PromisePegasus/worawich/Referense/hg38/gff3/Homo_sapiens.GRCh38.87.chr.gff3";
//        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/hg38_sim_3/";
//        String saveFilename = "hg38_FullNewMethod_Sim_alignmentResult_VariantReport";
        
//        int readLength = 24;
        int overlap = 4;
        byte percentMatch = 90;
        int coverageThreshold = 2;       
//        String saveFilenameCov = filename + "_VariantCoverageReport_match" + percentMatch;
        
        VariationResult varRes = SequenceUtil.analysisResultFromFileV3(fullPathSaveSortFile,merLength,overlap,percentMatch);
        varRes.createVariantReport();
 
        varRes.analyzeCoverageFusionV2();
        varRes.sortCoverageFusion();
//        varRes.writeVariantCoverageVirtualizeWithAnnotationReportToFile(filename2, gffFile, coverageThreshold, 'F');
        varRes.writeVariantSortedCoverageVirtualizeReportToFile(filename, coverageThreshold, 'F');
        varRes.writeVariantSortedCoverageReportToFile(filename, coverageThreshold, 'F',true);
        
        varRes.analyzeCoverageIndelV2();
        varRes.sortCoverageIndel();
//        varRes.writeVariantCoverageVirtualizeWithAnnotationReportToFile(filename2, gffFile, coverageThreshold, 'I');
        varRes.writeVariantSortedCoverageVirtualizeReportToFile(filename, coverageThreshold, 'I');
        varRes.writeVariantSortedCoverageReportToFile(filename, coverageThreshold, 'I',true);

        System.gc();
        /*****************************************************************************************/
    }
  
}
 