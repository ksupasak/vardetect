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
        String refPath = "/Users/worawich/Reference/hg38/hg38_filter_reIndex.fa";
        int numMer = 16;
        int numRead = 20;
        int readLen = 100;
        int readCoverage = 30;
        int diffL = 10000;
        int indelSize = 10;
        char variantType = 'T';                 // 4 variant type F=fusion L=large deeltion I=small insert D=small delete A=Large Insertion T=tandem Duplication
        char insertSNPFlag = 'F';               // T = true(insert) and F = false(not insert)
//        int minIndelSize = 10;
//        int maxIndelSize = 100;
        ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer); //runFile reference
        
        int maxMultiplier = 5;              // use to be the exponent of indelSizeBase
        int indelSizeBase = 10; 
        
//        InputSequence tempInSS = new InputSequence();
////        tempInSS = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomSingleType(ref, numRead, readLen, readCoverage, diffL, indelSize,filename,variantType);
//        tempInSS = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomSingleTypeFixRange(ref, numRead, readLen, readCoverage, minIndelSize, maxIndelSize, diffL, indelSize,filename,variantType);
//        System.out.println("done");
        
        for(int indelSizeMultiplier=3; indelSizeMultiplier<=maxMultiplier;indelSizeMultiplier++){        // indelSizeMultiplier start at 1
            long startTime = System.currentTimeMillis();
            
            int minIndelSize = (int)Math.pow(indelSizeBase, indelSizeMultiplier);                       // minIndelSize is indelSizeBase power by indelSizeMultiplier
            int maxIndelSize = (int)Math.pow(indelSizeBase, indelSizeMultiplier+1);                     // maxIndelSize is indelSizeBase power by indelSizeMultiplier+1
            
            
            String sampleFilename = "_tandemSim_"+variantType+"_SNP_"+insertSNPFlag+"_"+indelSizeMultiplier; 
//            ReferenceSequence ref = SequenceUtil.getReferenceSequence(refPath,numMer); //runFile reference
            /**
             * Simulate data
             */
            InputSequence input = new InputSequence();
    //        tempInSS = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomSingleType(ref, numRead, readLen, readCoverage, diffL, indelSize,filename,variantType);
            input = SimulatorUtil_WholeGene.simulateComplexWholeGeneRandomSingleTypeFixRange(ref, numRead, readLen, readCoverage, minIndelSize, maxIndelSize, diffL, indelSize,sampleFilename,variantType,insertSNPFlag);
            System.out.println("done");
            
            if(variantType == 'F'){
                break;
            }
            /*********************/
        }

    }
    
}
