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
public class TestIntersectRead {
    
    public static void main(String args[]) throws IOException{
        
//        String fasta1 = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/test_Extract_ABC/ERR718259_full_newRef_Mapped.fa";
//        String fasta2 = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/test_Extract_ABC/ERR718259_Mapped.fa";
//        
//        FastaUtil.createIntesectFastaFile(fasta1, fasta2);
        
        String fq1 = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/test_Extract_ABC/Test_extend_300/ERR718259_full_newRef_Mapped_unIntersect_ERR718259_unMapped_newRef.fq";
        String fq2 = "/Volumes/PromisePegasus/worawich/Download_dataset/TB_ERR718259/test_Extract_ABC/Test_extend_300/ERR718259_full_OriRef_Mapped.fq";
        
        FastqUtil.createIntesectFastqFile(fq1, fq2);
        
    }
}
