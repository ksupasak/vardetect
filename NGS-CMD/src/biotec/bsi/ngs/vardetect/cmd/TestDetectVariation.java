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
        
        String filename = "hg38_CHR16NewMethod_tha5.2_alignmentResult_Sorted";
        String path = "/Volumes/PromisePegasus/worawich/VMdev/dataScieneToolBox/projects/NGS/Thalassemia/";;
        String saveFilename = "hg38_simData_mul_alignmentResult_th5_VariantReport";
        String saveFilenameCov = "hg38_CHR16NewMethod_tha5.2_alignmentResult_VariantCoverageReport";
        

        VariationResult varRes = SequenceUtil.analysisResultFromFile(path+filename+".txt",18,100,4);
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
