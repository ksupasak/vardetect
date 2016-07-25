/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.inverseSequence;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 *
 * @author worawich
 */
public class NGSCMD5 {
    
    public static void main(String[] args) throws IOException { 
        long mask36Bit = 68719476735L;
        long mask2 = 268435455;
        String testBase = "101100001";
//        CharSequence testSeq = testBase;
//        String aon = "";
//        
//        for (int i = 0;i<testSeq.length();i++){
//            aon = aon + String.valueOf(testSeq.charAt((testSeq.length()-1)-i));
//        }
        //String invSeq = SequenceUtil.inverseSequence(testBase);
        //String outBase = SequenceUtil.createComplimentV2(invSeq);
        String revBin = new StringBuilder(testBase).reverse().toString();
        long revNum = Long.parseLong(revBin,2);
        long dummyNewMer = (~revNum)&mask36Bit;
        //System.out.println("Check inseq :\t" + testBase);
        //System.out.println("Check invseq :\t" + invSeq);
        //System.out.println("Check outseq :\t" + outBase);  
        System.out.println("Check outseq :\t" + testBase); 
        System.out.println("Check outseq :\t" + revBin);
        System.out.println("Check outseq :\t" + revNum);
        System.out.println("Check outseq :\t" + Long.toBinaryString(dummyNewMer));
        System.out.println("Check outseq :\t" + dummyNewMer);
        
        File f = new File("/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/hg19/chr21.bin"); //File object
        File compF = new File("/Users/worawich/VMdev/dataScieneToolBox/projects/NGS/hg19/chr21_comp.bin");
        DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
        DataInputStream compis = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
        int size = is.readInt();
        int compsize = compis.readInt();
        long[] list = new long[size];
        long[] complist = new long[compsize];
        for(int i=0;i<size;i++){
                  
            list[i] = is.readLong();
            
            
            long dummyMer = list[i]>>28;
            long dummyPos = list[i]&mask2;            
            
            for(int j=0;j<size;j++){
                complist[j] = compis.readLong();
                long compdummyMer = complist[j]>>28;
                long compdummyPos = complist[j]&mask2;
                if (dummyPos == compdummyPos){
                    System.out.println("Check dummyMer : " + dummyMer);
                    System.out.println("Check compdummyMer : " + dummyMer);
                    System.out.println("same position : " + compdummyPos + " and : " + compdummyPos);
                    System.out.println("Binary Check dummyMer : " + Long.toBinaryString(dummyMer));
                    System.out.println("Binary Check compdummyMer : " + Long.toBinaryString(compdummyMer));
                    break;
                }
            }
        }
        
    }
    
}
