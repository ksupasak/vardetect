/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;
import static biotec.bsi.ngs.vardetect.core.util.SVPUtility.modifyFile;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class RunModifyFile {
    
    public static void main(String[] args) throws IOException{
        
        String fileName = "/Volumes/PromisePegasus/worawich/Download_dataset/FASTQ_Taane_Martin/List_RD239_sample.txt";
        modifyFile(fileName);
    }
    
}
