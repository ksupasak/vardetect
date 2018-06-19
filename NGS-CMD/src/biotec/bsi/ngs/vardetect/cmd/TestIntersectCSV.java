/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.SVPUtility;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestIntersectCSV {
    
    public static void main(String[] args) throws IOException{
        
        String normalSample = "436";
        String tumorSample = "436T";
        
        String fileA = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+normalSample+"_sorted/Deletion/"+normalSample+"_sorted_deletion.csv";
        String fileB = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+tumorSample+"_sorted/Deletion/"+tumorSample+"_sorted_deletion.csv";
        String outPath = "/Users/worawich/Download_dataset/Ratina_cancer/summaryRes";
        String svType = "DEL";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+normalSample+"_sorted/Tandem/"+normalSample+"_sorted_tandem.csv";
        fileB = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+tumorSample+"_sorted/Tandem/"+tumorSample+"_sorted_tandem.csv";
        outPath = "/Users/worawich/Download_dataset/Ratina_cancer/summaryRes";
        svType = "TAN";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+normalSample+"_sorted/Interinsertion/"+normalSample+"_sorted_interinsertion.csv";
        fileB = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+tumorSample+"_sorted/Interinsertion/"+tumorSample+"_sorted_interinsertion.csv";
        outPath = "/Users/worawich/Download_dataset/Ratina_cancer/summaryRes";
        svType = "INTER";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+normalSample+"_sorted/Intrainsertion/"+normalSample+"_sorted_intrainsertion.csv";
        fileB = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+tumorSample+"_sorted/Intrainsertion/"+tumorSample+"_sorted_intrainsertion.csv";
        outPath = "/Users/worawich/Download_dataset/Ratina_cancer/summaryRes";
        svType = "INTRA";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+normalSample+"_sorted/Chimeric/"+normalSample+"_sorted_chimeric.csv";
        fileB = "/Users/worawich/Download_dataset/Ratina_cancer/res/"+tumorSample+"_sorted/Chimeric/"+tumorSample+"_sorted_chimeric.csv";
        outPath = "/Users/worawich/Download_dataset/Ratina_cancer/summaryRes";
        svType = "CHI";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
    }
    
}
