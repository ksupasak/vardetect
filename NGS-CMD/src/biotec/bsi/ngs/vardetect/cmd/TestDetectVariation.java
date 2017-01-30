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
        
        String filename = "hg38_simData_mul_alignmentResult_th5_PostLinuxSorted";       
        String path = "/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/test_sim_01/";
        String saveFilename = "hg38_simData_mul_alignmentResult_th5_VariantReport";
        String saveFilenameCov = "hg38_simData_mul_alignmentResult_th5_VariantCoverageReport";
        

        VariationResult varRes = SequenceUtil.analysisResultFromFile(path+filename+".txt",18,100,1);
        varRes.createVariantReport();
        varRes.analyzeCoverageFusion();
        varRes.writeVarianCoverageReportToFile(path, saveFilenameCov, 'F');
        varRes.writeVarianReportToFile(path, saveFilename);
        System.gc();
      
    }
  
}
