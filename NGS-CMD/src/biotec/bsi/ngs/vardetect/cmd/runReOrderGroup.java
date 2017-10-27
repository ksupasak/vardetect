/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.FastaUtil;
import biotec.bsi.ngs.vardetect.core.util.Report;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class runReOrderGroup {
    public static void main(String[] args) throws IOException{
        
        String sumReportFile = "/Volumes/PromisePegasus/worawich/Download_dataset/FASTQ_Taane_Martin/sum_CoverageReport_Indel.txt";
        Report.reOrderGroupNumber(sumReportFile);
//            
    }
}
