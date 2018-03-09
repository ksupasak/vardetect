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
import java.security.NoSuchAlgorithmException;

/**
 *
 * @author worawich
 */
public class TestUnIntersectRead {
    
    public static void main(String args[]) throws IOException, NoSuchAlgorithmException{
        // implement md5 check for sequence similarity (with sequence md5 check it can get rid of read that has the same sequence base but diferent in read name also)
        // fq1 should be the smallest file that contain read that we want to filter out
        
//        String fq1 = "/Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/NA12878-L1_S1_C.fq.gz";
//        String fq2 = "/Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/NA12878-L1_S1_map_preA.fq.gz";
//        String saveFileName = "/Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/testmd5unintersect.fq.gz";
        String fq1 = args[0];
        String fq2 = args[1];
        String saveFileName = args[2];
        
        FastqUtil.createUnIntesectFastqFile(fq1, fq2, saveFileName);
    }
    
}
