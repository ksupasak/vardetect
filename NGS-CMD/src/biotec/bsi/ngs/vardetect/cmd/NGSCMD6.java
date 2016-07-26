/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import java.io.IOException;

/**
 *
 * @author worawich
 */
public class NGSCMD6 {
    
    public static void main(String[] args) throws IOException { 
        long mask36Bit = 68719476735L;
        long mask2 = 268435455;
        long dummyMerPos = -9223372036854708401L;
//            System.out.println("Check fullcode dummyMerPos: " + dummyMerPos);
        long dummyMer = (dummyMerPos>>28)&mask36Bit;
        long dummyPos = dummyMerPos&mask2;

        /* Reconstruct compliment DNA sequence of whole chromosome */
//            System.out.println("Check mer sequemce befor compliment : " + dummyMer);
//            long dummyNewMer = (~dummyMer)&mask36Bit;
        System.out.println("Check binaryCode (input) : " + Long.toBinaryString(dummyMerPos));
        System.out.println("Check binaryCode (dummyMer) : " + Long.toBinaryString(dummyMer));
        System.out.println("Check binaryCode (dummyPos) : " + Long.toBinaryString(dummyPos));
        String binaryMer = Long.toBinaryString(dummyMer);
        int kmer = binaryMer.length()/2;
        System.out.println("Check binaryMer (before compliment) : " + binaryMer);
        System.out.println("Check dummyPos : " + dummyPos);
        
        
    }
    
}
