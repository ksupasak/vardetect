/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.FastaUtil;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class CreateSampleFromAlignResult {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here

        String resultFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_300517/dm6_300517_OP3_R1_alignResult_Sorted.txt";
        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/TCGA-75-5147-01A.unmapped.fa";


        FastaUtil.createSampleFromAlignResult(resultFile, sampleFile, 'c');

    }
}
