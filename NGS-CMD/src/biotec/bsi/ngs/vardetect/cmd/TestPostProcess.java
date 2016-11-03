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
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestPostProcess {
       public static void main(String[] args) throws IOException {
        // TODO code application logic here
        
        String savefilename = "test_New_structure_with_iniIndex_check";
       
    
        AlignmentResultRead readAlign = SequenceUtil.readAlignmentReportV2("/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/hg19"+savefilename+".txt");
        System.out.println("Begin create color array");
        Clustering.createColorArray(readAlign, 100, 18);
        System.out.println("Done create color array");
    }
  
}
