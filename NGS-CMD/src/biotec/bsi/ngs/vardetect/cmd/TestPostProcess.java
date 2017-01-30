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
public class TestPostProcess {
       public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
        String filename = "hg38_full_3661_mul_alignmentResult_th5_part";       
        String path = "/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/test_sim_01/";
        String saveFileName = "hg38_3661_alignmentResult_forLinuxSort";
        String filename3661 = "hg38_simData_mul_alignmentResult_th5_part";
        String filenameERR = "hg38_simData_mul_alignmentResult_th5_ERR_part";
        String saveFileNameERR = "hg38_simData_mul_alignmentResult_th5_ERR_forLinuxSortV2";
        String saveFileName3661 = "hg38_simData_mul_alignmentResult_th5_forLinuxSort";
        String saveFileNameForPostProcess = "hg38_simData_mul_alignmentResult_th5_PostLinuxSortedV2";
        int numPart = 1;
        
        for(int i=1;i<=numPart;i++){
            
            AlignmentResultRead readAlign = SequenceUtil.readAlignmentReportV2(path+filename3661+i+".txt",100,18);
            System.out.println("Begin create color array");
            Clustering.createColorArray(readAlign, 100, 18);        
            System.out.println("Done create color array");
            readAlign.writeSortedCutColorResultToPathInFormatForLinuxSort(path, saveFileName3661, "txt");
//            VariationResult varRes = SequenceUtil.analysisResultFromFile(path+saveFileNameForPostProcess+".txt",18,100);
            readAlign = null;
            System.gc();
       }
    }
  
}
