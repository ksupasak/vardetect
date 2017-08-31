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
public class AnalyzeFusionFromIntersectLASTResult {
    public static void main(String args[]) throws IOException{
        
         /**
          * Plug in intersect result File and threshold for filter low coverage  
          */
         
        String intersectLASTResult1 = args[0];
        int threshold = Integer.parseInt(args[1]);

//        String intersectLASTResult1 = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/LAST_result_eclipse/unmapped_cancer/TCGA-55-6543-o1A_unmapped_NC_001414_hg38Filter_LASTIntesectResult.txt";        
//        int threshold = 1;        
        LASTUtils.analyzeFusionFromIntersecResult(intersectLASTResult1, threshold);
                
    }
}
