/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.VariationResult;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestPostProcessV2 {
    
    public static void main(String[] args) throws IOException{
        
        String inputFile = "/Users/worawich/Download_dataset/Ratina_cancer/277_sorted_unmap.filter.out";
        String saveFile = "/Users/worawich/Download_dataset/Ratina_cancer/277_sorted_unmap.filter";
        VariationResult varRes = SequenceUtil.readVersion2AlignmentResult(inputFile);
        varRes.analyzeCoverage();
        varRes.classifyRoughSVType();
        varRes.writeStructureVariantV2SortedCoverageReportToFile(saveFile, 5);
        varRes.writeStructureVariantV2SortedCoverageGroupInfoReportToFile(saveFile, 5);
    }
}