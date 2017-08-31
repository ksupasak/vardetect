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
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer); //runFile hg19.fa
        
        InputSequence tempInSS = new InputSequence();
        tempInSS = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomMixed(ref, 200, 100, 10, 10000, 9,"hg38_sim_2");
        System.out.println("done");

    }
    
}
