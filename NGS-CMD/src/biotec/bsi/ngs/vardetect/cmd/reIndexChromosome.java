/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.FastaUtil;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;

/**
 *
 * @author worawich
 */
public class reIndexChromosome {
    
    public static void main(String[] args) throws IOException{
        
        String fastaFile = "/Volumes/PromisePegasus/worawich/Referense/TB_reference/H37Rv_NC_000962.3.txt";

        FastaUtil.reIndexChrNameFastaFile(fastaFile);
//            
    }
    
}
