/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.FastaUtil;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class AnalyzeNonVariant {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        /**
         * This function will analyze the coverage of each pattern exist in each read
         * input file would be in format alignResult_Sort.txt
         */
        
        String nonVariantFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_050417/OP3_S3_dm6_miRNA_mer12_alignResult_forLinuxSort.txt";
        
        int merLength = 12;
        int coverageThreshold = 1;
        SequenceUtil.analysisNonVariantResultFromFile(nonVariantFile, merLength, coverageThreshold);
 
    }
}
