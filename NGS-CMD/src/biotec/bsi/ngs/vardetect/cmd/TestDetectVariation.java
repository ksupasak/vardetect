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
        
        String filename = "/Volumes/PromisePegasus/worawich/Download_dataset/Thalasemia_Data/3661/3661_unmapped_alignResult_Sort.txt";
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

        varRes.analyzeCoverageFusion();
        varRes.writeVariantCoverageReportToFile(filename, coverageThreshold, 'F');
        varRes.writeVariantCoverageVirtualizeReportToFile(filename, coverageThreshold, 'F');
        varRes.analyzeCoverageIndel();
        varRes.writeVariantCoverageReportToFile(filename, coverageThreshold, 'I');
        varRes.writeVariantCoverageVirtualizeReportToFile(filename, coverageThreshold, 'I');

        System.gc();
      
    }
  
}
