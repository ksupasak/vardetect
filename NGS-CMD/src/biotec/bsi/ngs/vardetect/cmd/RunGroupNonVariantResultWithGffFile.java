/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class RunGroupNonVariantResultWithGffFile {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        /**
         * This MD will run grouping the alignment result 
         * it require result that passing post process and already sorted
         * It group the result with gff3 file or group it by annotation
         */
        
        String resultFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_300517/dm6_300517_OP4_R2_alignResult_Sorted.txt";
        String gffFile = "/Volumes/PromisePegasus/worawich/Referense/drosophila/ensemble/Drosophila_melanogaster.BDGP6.89.gff3_2_cutGene.txt";
        
        SequenceUtil.groupNonVariantResultWithGffFile(gffFile, resultFile);
 
    }
}
