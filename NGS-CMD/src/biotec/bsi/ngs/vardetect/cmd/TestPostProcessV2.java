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
        
//        String inputFile = "/Users/worawich/Download_dataset/verayuth/res/T102_dad.recal/T102_dad.recal_unmap_filter.out";
//        String inputFile = "/Users/worawich/Download_dataset/DupDel/res_rmDup/SEAdel/SEAdel_unmap_filter.out";
//        String inputFile = "/Users/worawich/Download_dataset/SLE/SLE_bam/res_rmDup/SLE-5.recal/del_del_filter/SLE-5.recal_unmap_filter.out";
        String inputFile = "/Users/worawich/Download_dataset/Ratina_cancer/277T_sorted_unmap_filter.out";
//        String saveFile = "/Users/worawich/Download_dataset/SLE/SLE_bam/res_rmDup/SLE-5.recal/SLE-5.recal_unmap_filter";
        String refFile = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa";
        String refFaIdx = "/Users/worawich/Reference/hg19_SVP2/hg19_main.fa.fai";
        String gffFile = "/Users/worawich/Reference/hg19_SVP2/Homo_sapiens.GRCh37.87.chr.gff3";
//        String bamFile = "/Users/worawich/Download_dataset/verayuth/T102_dad.recal.bam";
//        String bamFile = "/Users/worawich/Download_dataset/DupDel/res_rmDup/SEAdel/SEAdel_removeDup.bam";
        String bamFile = "/Users/worawich/Download_dataset/Ratina_cancer/277T_sorted.bam";
//        String bamFile = "/Users/worawich/Download_dataset/SLE/SLE_bam/SLE-5.recal.bam";
        String samtoolsDirectory = "/usr/local/bin/samtools";
        
        VariationResult varRes = SequenceUtil.readVersion2AlignmentResult(inputFile);
        varRes.createRefIndex(refFaIdx);
        varRes.analyzeCoverage();
        
        int numPerGroup = 3;
        int maxExtend = 100;
        int numMer = 16;
//        varRes.classifyRoughSVType();
//        varRes.writeStructureVariantV2SortedCoverageReportToFile(saveFile, 5);
//        varRes.writeStructureVariantV2SortedCoverageGroupInfoReportToFile(saveFile, 5);
//        varRes.writeStructureVariantV2SortedCoverageReportWithAnnotationToFile(saveFile, gffFile, refFaIdx, 5, 16);
//        varRes.writeStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationToFile(saveFile, gffFile, refFaIdx, 5, 16);
        
        varRes.classifyPreciseSVType(numPerGroup);
        varRes.identifyCorrectness(bamFile,samtoolsDirectory);
        varRes.sortInsertionByInsertJunctiontLowtoHigh();
        varRes.FilterIntraInsertion();  // Del Del filter
        
        varRes.writePreciseStructureVariantV2SortedCoverageGroupInfoReportWithAnnotationExcel(inputFile, gffFile, refFaIdx, numPerGroup, numMer);
        varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refFile, refFaIdx, gffFile, "TD", numPerGroup, maxExtend, numMer);
        varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refFile, refFaIdx, gffFile, "D", numPerGroup, maxExtend, numMer);
        varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refFile, refFaIdx, gffFile, "IA", numPerGroup, maxExtend, numMer);
        varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refFile, refFaIdx, gffFile, "IE", numPerGroup, maxExtend, numMer);
        varRes.writeVisualizePreciseSVTypeWithAnnotation(inputFile, refFile, refFaIdx, gffFile, "CH", numPerGroup, maxExtend, numMer);
//        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "TD", 5, 100);
//        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "D", 5, 100);
//        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "IA", 5, 100);
//        varRes.writeVisualizePreciseSVType(inputFile, refFile, refFaIdx, "IE", 5, 100);
//        varRes.writeStructureVariantV2EuclideanDistanceTable(saveFile, 5);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "TD", 5, 400);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "ID", 5, 400);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "IC", 5, 400);
//        varRes.createReferenceFromNovelIndelResult_VariationV2(inputFile, refFile, refFaIdx, "IT", 5, 400);
    }
}
