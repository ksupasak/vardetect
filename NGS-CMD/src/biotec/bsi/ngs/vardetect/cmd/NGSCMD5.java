/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.inverseSequence;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class NGSCMD5 {
    
    public static void main(String[] args) throws IOException { 
        
        String testBase = "CAAAATGCTTCGTAAGTCAAATTTGGGCAGTG";
//        CharSequence testSeq = testBase;
//        String aon = "";
//        
//        for (int i = 0;i<testSeq.length();i++){
//            aon = aon + String.valueOf(testSeq.charAt((testSeq.length()-1)-i));
//        }
        String invSeq = SequenceUtil.inverseSequence(testBase);
        String outBase = SequenceUtil.createCompliment(invSeq);
        
        System.out.println("Check inseq :\t" + testBase);
        System.out.println("Check invseq :\t" + invSeq);
        System.out.println("Check outseq :\t" + outBase);  
        
    }
    
}
