/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.FastaUtil;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class RunCreateSampleFromCoverageReport {
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
//        String sampleFile = "/Volumes/PromisePegasus/worawich/Download_dataset/FASTQ_Taane_Martin/ERR846998_2.fastq";
//        String coverageFile = "/Volumes/PromisePegasus/worawich/Download_dataset/FASTQ_Taane_Martin/ERR846998_2_H37RV_alnRes_part_match70_CoverageReport_Indel.txt";
        String coverageReportFileName = args[0];
        String sampleFile = args[1];

        FastaUtil.createSampleFromVariantCoverageReport(coverageReportFileName, sampleFile);
    }
}
