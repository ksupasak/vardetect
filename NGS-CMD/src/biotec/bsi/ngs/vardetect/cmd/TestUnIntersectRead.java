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
        
//        String fq1 = "/Users/worawich/Download_dataset/TB/SampleForTestScript/debug/ERR718193_C.fq";
//        String fq2 = "/Users/worawich/Download_dataset/TB/SampleForTestScript/debug/ERR718193_map_preA.fq.gz";
//        String saveFileName = "/Users/worawich/Download_dataset/TB/SampleForTestScript/debug/testmd5unintersect.fq";
        String fq1 = args[0];
        String fq2 = args[1];
        String saveFileName = args[2];
        
        FastqUtil.createUnIntesectFastqFile(fq1, fq2, saveFileName);
    }
    
}
