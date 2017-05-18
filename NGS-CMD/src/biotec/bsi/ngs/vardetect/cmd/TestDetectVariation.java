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

/**
 *
 * @author worawich
 */
public class TestDetectVariation {
       public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
        String filename = "hg38_sim_4thread_th5_alignmentResult_LongRead_Sorted";
        String path = "/Volumes/PromisePegasus/worawich/Download_dataset/SimulateData/test_sim_01/";
//        String saveFilename = "hg38_FullNewMethod_Sim_alignmentResult_VariantReport";
        String saveFilenameCov = "hg38_sim_4thread_th5_alignmentResult_LongRead_VariantCoverageReport";
//        int readLength = 24;
        int merLength = 18;
        int overlap = 4;

        VariationResult varRes = SequenceUtil.analysisResultFromFileV3(path+filename+".txt",merLength,overlap);
        varRes.createVariantReport();
        varRes.analyzeCoverageFusion();
        varRes.writeVariantCoverageReportToFile(path, saveFilenameCov, 'F');
        varRes.analyzeCoverageIndel();
        varRes.writeVariantCoverageReportToFile(path, saveFilenameCov, 'I');
        varRes.analyzeCoverageSNP();
        varRes.writeVariantCoverageReportToFile(path, saveFilenameCov, 'S');
//        varRes.writeVarianReportToFile(path, saveFilename);
        System.gc();
      
    }
  
}
