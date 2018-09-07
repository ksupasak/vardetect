/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.CombineReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class runPreIndexReference {
    
    public static void main(String[] args) throws IOException, FileNotFoundException, InterruptedException{
        String refPath = args[0];
        int numMer = Integer.parseInt(args[1]);
        int numThread = Integer.parseInt(args[2]);
        
        CombineReferenceSequence ref = SequenceUtil.getCombineReferenceSequence(refPath,numMer,numThread); //runFile hg19.fa
        ref.setMaximumDuplicatePattern(10);
        ref.prepare();
    }
    
}
