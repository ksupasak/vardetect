/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_WholeGene;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class RunSimulateSample {
    
    public static void main(String args[]) throws IOException{
        String refPath = "/Volumes/PromisePegasus/worawich/Referense/hg38/hg38_filter.fa";
        int numMer = 18;
        int numRead = 100;
        int readLen = 100;
        int readCoverage = 30;
        int diffL = 10000;
        int indelSize = 10;
        String filename = "hg38_sim_L";
        char variantType = 'L';                 // 4 variant type F=fusion L=large indel I=insert D=delete 
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer); //runFile reference
        
        InputSequence tempInSS = new InputSequence();
        tempInSS = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomSingleType(ref, numRead, readLen, readCoverage, diffL, indelSize,filename,variantType);
        System.out.println("done");

    }
    
}
