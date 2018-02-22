/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.SVPUtility;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;

/**
 *
 * @author worawich
 */
public class testModifyBam {
    public static void main(String[] args) throws NoSuchAlgorithmException, IOException{
        String mainBamFile = "/Users/worawich/Download_dataset/TB/SampleForTestScript/debug/ERR718193.bam";
        String targetBamFile = "/Users/worawich/Download_dataset/TB/SampleForTestScript/debug/ERR718193_C_map.bam";
        String outputFormat = "sam";
//        
//        SVPUtility.addUnmapFlagToBam(mainBamFile, targetBamFile, outputFormat);
        
//        String mainBamFile = args[0];
//        String targetBamFile = args[1];
//        String outputFormat = args[2];
        
        SVPUtility.addUnmapFlagToBam(mainBamFile, targetBamFile, outputFormat);
                
    }
}
