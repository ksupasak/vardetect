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
public class TestContig {
    public static void main(String[] args) throws IOException{
    
        String[] filename = {"/Volumes/PromisePegasus/worawich/Download_dataset/crab_transscriptome/referencecrab/fromgigadb/crab"};
        
        
        FastaUtil.createReferenceFromContig(filename);
    }
}
