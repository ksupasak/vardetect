/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.VariationResult;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestPostProcessV2 {
    
    public static void main(String[] args) throws IOException{
        
//        String inputFile = "/Users/worawich/Download_dataset/SLE/SLE_bam/Bam_dedup/test_noIndel/res_testnewSort/SLE-9_rmdup/SLE-9_rmdup_unmap_filter.out";
//        String saveFile = "/Users/worawich/Download_dataset/SLE/SLE_bam/Bam_dedup/test_noIndel/res_testnewSort/SLE-9_rmdup/SLE-9_rmdup_unmap_filter";
//        String refFile = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa";
//        String refFaIdx = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa.fai";
//        String gffFile = "/Users/worawich/Reference/hg19_SVP2/Homo_sapiens.GRCh37.87.chr.gff3";
        
//        String inputFile = "/Users/worawich/Download_dataset/Ratina_cancer/432T_sorted_unmap_filter.out";
//        String saveFile = "/Users/worawich/Download_dataset/Ratina_cancer/432T_sorted_unmap_filter";
//        String refFile = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa";
//        String refFaIdx = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa.fai";
//        String gffFile = "/Users/worawich/Reference/hg19_SVP2/Homo_sapiens.GRCh37.87.chr.gff3";
        

        String inputFile = "/Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/Result_firstVersion_nonSeparateJunction/NA12878-L1_S1_unmap_filter.out";
        String saveFile = "/Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/Result_firstVersion_nonSeparateJunction/NA12878-L1_S1_unmap_filter";
        String refFile = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa";
        String refFaIdx = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa.fai";
        String gffFile = "/Users/worawich/Reference/hg19_SVP2/Homo_sapiens.GRCh37.87.chr.gff3";
        
        String bamFile = "/Users/worawich/Download_dataset/HiSeq_X_Ten_TruSeq_Nano_4_replicates_of_NA12878_-8998991/BWA_Whole_Genome_Sequencing_NA12878_L1-10532644/NA12878_L1-ds_d51e1a3660b647f09677f58400c86764/rmDup_Bam/NA12878-L1.rmdup.bam";
        String samtoolsDirectory = "/usr/local/bin/samtools";
        
        VariationResult varRes = SequenceUtil.readVersion2AlignmentResult(inputFile);
        varRes.createRefIndex(refFaIdx);
        varRes.analyzeCoverage();
//        varRes.classifyRoughSVType();
//        varRes.writeStructureVariantV2SortedCoverageReportToFile(saveFile, 5);
//        varRes.writeStructureVariantV2SortedCoverageGroupInfoReportToFile(saveFile, 5);
//        varRes.writeStructureVariantV2SortedCoverageReportWithAnnotationToFile(saveFile, gffFile, refFaIdx, 5, 16);
//        varRes.writeStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationToFile(saveFile, gffFile, refFaIdx, 5, 16);
        
        varRes.classifyPreciseSVType(5);
        //varRes.identifyCorrectness(bamFile,samtoolsDirectory);
        //varRes.writePreciseStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationToFile(saveFile, gffFile, refFaIdx, 5, 16);
        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "TD", 5, 100);
        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "D", 5, 100);
        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "IA", 5, 100);
        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "IE", 5, 100);
//        varRes.writeStructureVariantV2EuclideanDistanceTable(saveFile, 5);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "TD", 5, 400);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "ID", 5, 400);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "IC", 5, 400);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "IT", 5, 400);
    }
}
