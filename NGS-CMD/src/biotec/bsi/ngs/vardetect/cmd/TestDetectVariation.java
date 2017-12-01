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
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author worawich
 */
public class TestDetectVariation {
       public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
//        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_Error/TB_Error_RD239_Sort.txt";
        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/Test_ERR718259_genRef/ERR718259_unMap_Sort.txt";
//        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/Error/error_cancer.txt";
        String gffFile = "/Volumes/PromisePegasus/worawich/Referense/hg38/gff3/Homo_sapiens.GRCh38.87.chr.gff3";
//        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/hg38_sim_3/";
//        String saveFilename = "hg38_FullNewMethod_Sim_alignmentResult_VariantReport";
        
//        int readLength = 24;
        int merLength = 18;
        int overlap = 4;
        byte percentMatch = 90;
        int coverageThreshold = 2;       
//        String saveFilenameCov = filename + "_VariantCoverageReport_match" + percentMatch;
        
        VariationResult varRes = SequenceUtil.analysisResultFromFileV3(filename,merLength,overlap,percentMatch);
        varRes.createVariantReport();
//        varRes.analyzeCoverageFusion();
//        varRes.writeVariantCoverageReportToFile(path, saveFilenameCov, coverageThreshold, 'F');
//        varRes.analyzeCoverageIndel();
//        varRes.writeVariantCoverageReportToFile(path, saveFilenameCov, coverageThreshold, 'I');
//        varRes.analyzeCoverageSNP();
//        varRes.writeVariantCoverageReportToFile(path, saveFilenameCov, coverageThreshold, 'S');
//        varRes.writeVarianReportToFile(path, saveFilename);

//        varRes.analyzeCoverageFusion();
//        varRes.writeVariantCoverageReportToFile(filename, coverageThreshold, 'F');
//        varRes.writeVariantCoverageVirtualizeReportToFile(filename, coverageThreshold, 'F');
//        varRes.analyzeCoverageIndel();
//        varRes.writeVariantCoverageReportToFile(filename, coverageThreshold, 'I');
//        varRes.writeVariantCoverageVirtualizeReportToFile(filename, coverageThreshold, 'I');
        
        varRes.analyzeCoverageFusionV2();
        varRes.sortCoverageFusion();
//        varRes.writeVariantCoverageVirtualizeWithAnnotationReportToFile(filename, gffFile, coverageThreshold, 'F');
        varRes.writeVariantSortedCoverageVirtualizeReportToFile(filename, coverageThreshold, 'F');
        varRes.writeVariantSortedCoverageReportToFile(filename, coverageThreshold, 'F',true);             // last argunement set true tell us it generate read name only report and normal indel coverage report
        varRes.analyzeCoverageIndelV2();
        varRes.sortCoverageIndel();
//        varRes.writeVariantCoverageVirtualizeWithAnnotationReportToFile(filename, gffFile, coverageThreshold, 'I');
        varRes.writeVariantSortedCoverageVirtualizeReportToFile(filename, coverageThreshold, 'I');
        varRes.writeVariantSortedCoverageReportToFile(filename, coverageThreshold, 'I',true);         
        
        String refFaIndex = "/Volumes/PromisePegasus/worawich/Referense/TB_reference/H37Rv_NC_000962_reIndex.fa.fai";
        String refFile = "/Volumes/PromisePegasus/worawich/Referense/TB_reference/H37Rv_NC_000962_reIndex.fa";
        String sampleFaIndex = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/ERR718259_unmaped.fa.fai";
        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/ERR718259_unmaped.fa";
        varRes.createReferenceFromNovelIndelResult(filename, refFile, refFaIndex, sampleFile, sampleFaIndex, "LI", 10, 100);
        // Export new indel type into bed 12 format for create reference (SI = small insertion, SD = small deletion, LI=large insertion)
//        varRes.exportNewIndelEventToBed12(filename, "LI", 10);
        System.gc();
      
    }
  
}
