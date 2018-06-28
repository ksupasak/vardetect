/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.SVPUtility;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestIntersectCSV {
    
    public static void main(String[] args) throws IOException{
        
        String normalSample = "277";
        String tumorSample = "277T";
        String parentPath = "/Users/worawich/Download_dataset/Ratina_cancer/res_with_p54_del_del_filter";
        String savePath = "/Users/worawich/Download_dataset/Ratina_cancer/res_intersect"+File.separator+normalSample+"_"+tumorSample;
        
        File file = new File(savePath);
        
        if(!file.getParentFile().exists()){
            try{
                file.getParentFile().mkdirs();
            } 
            catch(SecurityException se){
                System.out.println("Can not create directory : " + file.getParentFile());
            }
        }
        
        if(!file.getAbsoluteFile().exists()){
            try{
                file.getAbsoluteFile().mkdirs();
            } 
            catch(SecurityException se){
                System.out.println("Can not create directory : " + file.getAbsoluteFile());
            }
        }
        
        
        
        
        String fileA = parentPath+File.separator+normalSample+"_sorted/Deletion/"+normalSample+"_sorted_deletion.csv";
        String fileB = parentPath+File.separator+tumorSample+"_sorted/Deletion/"+tumorSample+"_sorted_deletion.csv";
        String outPath = savePath;
        String svType = "DEL";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = parentPath+File.separator+normalSample+"_sorted/Tandem/"+normalSample+"_sorted_tandem.csv";
        fileB = parentPath+File.separator+tumorSample+"_sorted/Tandem/"+tumorSample+"_sorted_tandem.csv";
        outPath = savePath;
        svType = "TAN";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = parentPath+File.separator+normalSample+"_sorted/Interinsertion/"+normalSample+"_sorted_interinsertion.csv";
        fileB = parentPath+File.separator+tumorSample+"_sorted/Interinsertion/"+tumorSample+"_sorted_interinsertion.csv";
        outPath = savePath;
        svType = "INTER";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = parentPath+File.separator+normalSample+"_sorted/Intrainsertion/"+normalSample+"_sorted_intrainsertion.csv";
        fileB = parentPath+File.separator+tumorSample+"_sorted/Intrainsertion/"+tumorSample+"_sorted_intrainsertion.csv";
        outPath = savePath;
        svType = "INTRA";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
        
        fileA = parentPath+File.separator+normalSample+"_sorted/Chimeric/"+normalSample+"_sorted_chimeric.csv";
        fileB = parentPath+File.separator+tumorSample+"_sorted/Chimeric/"+tumorSample+"_sorted_chimeric.csv";
        outPath = savePath;
        svType = "CHI";
        
        SVPUtility.fastIntersectCSVReport(fileA, fileB, outPath, svType);
    }
    
}
