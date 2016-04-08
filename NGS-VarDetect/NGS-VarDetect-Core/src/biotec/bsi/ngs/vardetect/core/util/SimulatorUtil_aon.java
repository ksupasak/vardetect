/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import java.util.Vector;

/**
 *
 * @author Worawich
 */
public class SimulatorUtil_aon {
    
    public static InputSequence simulateWholeGene(ReferenceSequence ref, int num_read, int ln_read){
        
        InputSequence is = new InputSequence();
        Vector<ChromosomeSequence> chrs = ref.getChromosomes();
        
        System.out.println(chrs);
        
        
        
        
        
        
        
        return is;
    }
    
}
