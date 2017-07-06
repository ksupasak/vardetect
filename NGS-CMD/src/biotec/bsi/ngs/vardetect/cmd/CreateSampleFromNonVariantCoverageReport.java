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
public class CreateSampleFromNonVariantCoverageReport {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/TCGA-75-5147-01A.unmapped.fa";
        String coverageFile = "/Volumes/PromisePegasus/worawich/Download_dataset/cancer/unmaped_cancer/TCGA_75_5147_virus_alignResult_Sort_nonVariantCoverage.txt";

        FastaUtil.createSampleFromNonVariantCoverageReport(coverageFile, sampleFile);
    }
}
