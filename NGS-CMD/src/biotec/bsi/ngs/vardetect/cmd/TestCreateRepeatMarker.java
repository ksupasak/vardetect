/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class TestCreateRepeatMarker {
    public static void main(String[] args) throws IOException {
        String refPath = "/Volumes/PromisePegasus/worawich/VMdev/dataScieneToolBox/projects/NGS/hg38/hg38_filter.fa";
        int numMer = 18;
        
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer);   
        
    }    
    
}
