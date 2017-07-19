/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.alignment.AlignerFactory;
import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResult;
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.ClusterGroup;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.Clustering;
import biotec.bsi.ngs.vardetect.core.util.FastaUtil;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_WholeGene;
import biotec.bsi.ngs.vardetect.core.util.VisualizeResult;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */

public class TestBedCMD {
    
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
      

        String fastaFile = "/Volumes/PromisePegasus/worawich/Referense/NC_001494.1_Jaagsiekte_sheep_retrovirus/sequence.fasta";
        String resultFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_300517/dm6_300517_OP4_R2_alignResult_Sorted.txt";
        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/Micro_RNA/NGS_result_300517/OP4_S4_L001_R2_001_cutadapt_filter.fa";
        String nonVariantFile = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/TCGA_75_5147_virus_alignResult_Sort.txt";
        String coverageFile = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/TCGA_75_5147_virus_alignResult_Sort_nonVariantCoverage.txt";
//          int startPoint = 0;
//          int length = 15;
//          SequenceUtil.truncateFastaFIles(fastaFile, startPoint, length);
//          FastaUtil.filterSampleFile(fastaFile);
//          SequenceUtil.miRNASeparator(fastaFile);
//            FastaUtil.reIndexChrNameFastaFile(fastaFile);
        FastaUtil.createSampleFromAlignResult(resultFile, sampleFile, 'c');
//        SequenceUtil.analysisNonVariantResultFromFile(nonVariantFile, 18, 2);
//        FastaUtil.createSampleFromNonVariantCoverageReport(coverageFile, sampleFile);
        
//        FastaUtil.separateContigToFile(sampleFile);
    }

}
