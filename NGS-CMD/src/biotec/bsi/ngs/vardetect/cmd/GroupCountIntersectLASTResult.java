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
public class GroupCountIntersectLASTResult {
     public static void main(String args[]) throws IOException{
        
         /**
          * recommend that the lastResult1 must have lesser sample than LASTResult2 
          */
         
        String LASTResult1 = args[0];
        String LASTResult2 = args[1];
        
//        String LASTResult1 = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/LAST_result_eclipse/unmapped_cancer/TCGA-55-6543-o1A.unmapped_NC_001414.maf";
//        String LASTResult2 = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/LAST_result_eclipse/unmapped_cancer/LAST_hg38_Result/TCGA-75-5147-01A.unmapped_hg38_db.maf";
        
        LASTUtils.countIntersectFromLASTResult(LASTResult1, LASTResult2);
                
    }
}
