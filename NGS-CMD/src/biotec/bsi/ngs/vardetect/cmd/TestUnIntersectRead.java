/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.FastaUtil;
import biotec.bsi.ngs.vardetect.core.util.FastqUtil;
import biotec.bsi.ngs.vardetect.core.util.SamUtil;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestUnIntersectRead {
    
    public static void main(String args[]) throws IOException{
        
        String fq1 = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/test_Extract_ABC/Test_extend_300/potentialB/ERR718259_unMapped_newRef_unIntersect_ERR718259_full_CIGARCUT_newRef_Mapped.fq";
        String fq2 = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/test_Extract_ABC/Test_extend_300/ERR718259_full_OriRef_Mapped.fq";
        
        FastqUtil.createUnIntesectFastqFile(fq1, fq2);
    }
    
}
