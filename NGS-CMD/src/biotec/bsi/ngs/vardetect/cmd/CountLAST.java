/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.LASTUtils;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class CountLAST {
    
    public static void main(String args[]) throws IOException{
        
        String LASTFileName = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/LAST_result_eclipse/unmapped_cancer/TCGA-97-7938-01A.unmapped_hg38_filter_db.maf";
        LASTUtils.countSampleFromLASTResult(LASTFileName);
                
    }
}
