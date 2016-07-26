/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.alignment;

import biotec.bsi.ngs.vardetect.core.Aligner;
import biotec.bsi.ngs.vardetect.core.AlignmentResult;
import biotec.bsi.ngs.vardetect.core.AlignmentResultRead;
import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.MerRead;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import static biotec.bsi.ngs.vardetect.core.util.SequenceUtil.encodeSerialChromosomeSequenceV3;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author soup
 */
public class BinaryAligner implements Aligner{

    ReferenceSequence ref;
    
    
    public void setReferenceSequence(ReferenceSequence ref){
        this.ref = ref;
        
    }
    

   
    public AlignmentResultRead align(ReferenceSequence ref, InputSequence input) {
        
        this.setReferenceSequence(ref);
        return align(input);
        
    }
    
    public AlignmentResultRead align(InputSequence input){
        
        
        AlignmentResult res = new AlignmentResult(input);
        AlignmentResultRead alinResult = new AlignmentResultRead();
        
        int mer = 18;
        
        
        Enumeration<ChromosomeSequence> chrs = ref.getChromosomes().elements();

        while(chrs.hasMoreElements()){
            try {
                
                
                ChromosomeSequence chr = chrs.nextElement();
                Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
                System.out.println("reading .. "+chr.getName()+"");

                EncodedSequence encoded = encodeSerialChromosomeSequenceV3(chr);
                while(seqs.hasMoreElements()){
                    ShortgunSequence seq = seqs.nextElement();
                    Map<Long,long[]> merMap = new HashMap();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();
                    
                    System.out.print(chr.getName()+" + strand\t"); 
                    
                    
                    for(int i=0;i<(s.length()-mer)+1;i++){
                        int index = i;
                        String sub = s.substring(i, i+mer);
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){
                            m = m<<28;
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2(m);
//                            long pos2[] = encoded.fullAlign(m);
                            long pos = -1;
                            if(pos2!=null&&pos2.length>0){
                                pos = pos2[0];
                                pos = pos2.length;
                                //merMap = res.addResult(m, chr.getChrNumber(), pos2);
                            }
                            
                            int idx = (int) (pos-i);
                            if(pos<0){
                              idx = 0;
                            }
//                            System.out.println(""+chr.getName()+" "+sub+" "+sub.length()+" : "+m+" pos : "+pos+" idx : "+idx);
                            
                            System.out.print("\t"+pos);
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("Before Mermap this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            merMap = res.addResult(m, chr.getChrNumber(), pos2); // Still confuse After pss this line all member in pos2 will chang from 28 bit of position to 36 bit of chr:Pos Hoe did it happen !!!
                            // But it work!
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("After mer map this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            res.addResultV2(m, chr.getChrNumber(), pos2, seq.getReadName());
                            //System.out.println("Check seq length" + seq.getShortgunLength());
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Stroe in object)---------------------------------------------*/
                            if (seq.getMerReadSize() < totalMer){
                                MerRead merRead = new MerRead();
                                
                                merRead.addMatchResultStrand(m, pos2, index, chr.getChrNumber());
                                seq.addMerRead(merRead);
//                                System.out.println("(First Time) Size Mer Read check: " + seq.getMerReadSize());
                            }else{
//                                System.out.println("(Other time) Size Mer Read check: " + seq.getMerReadSize());
                                if (index<seq.getMerReadSize()){
//                                    System.out.println("Check Index: " + index);
                                    MerRead merRead = seq.getMerRead().get(index); // call back existing merRead to update
                                    merRead.addMatchResultStrand(m, pos2, index, chr.getChrNumber());
                                    seq.addMerReadByIndex(index,merRead);
//                                    System.out.println("(Other time) Should be constant: Size Mer Read check: " + seq.getMerReadSize());
                                }                                   
                            }
                            /*-----------------------------------------------------------------------------------------------------------*/
                            /*************************************************************************************************************/
                            
                        }
                    //System.out.println(" This mer Map check: "+ (merMap == null));
                    res.createMap(seq.getReadName(), merMap);                    
                    }
                     System.out.println();
                    /* New Implement Part */
                    //seq.countAlignmentData(); // Create Alignment count data before change ShortgunSequence
                    
                    /*--------------------*/
                }
                
                /*-------------------- Do compliment alignment -------------------------------*/
                Enumeration<ShortgunSequence> seqsComp = input.getInputSequence().elements();
                while(seqsComp.hasMoreElements()){
                    ShortgunSequence seq = seqsComp.nextElement();
                    Map<Long,long[]> merMap = new HashMap();
                    
//                    System.out.println(""+chr.getName()+" "+encoded.getMers().length);
                    
                    String s = seq.getSequence();
                    String invSeq = SequenceUtil.inverseSequence(s);
                    String compSeq = SequenceUtil.createComplimentV2(invSeq);
                    System.out.println("******Input Sequence check " + compSeq);
                    System.out.print(chr.getName()+" - strand\t"); 
                    
                    
                    for(int i=0;i<(compSeq.length()-mer)+1;i++){
                        int index = i;
                        String sub = compSeq.substring(i, i+mer);
                        //System.out.println("check sub length"+sub.length());
                        long m = SequenceUtil.encodeMer(sub, mer);
//                        System.out.println(""+sub+" "+sub.length()+": "+m);
                        if(m!=-1){
                            m = m<<28;
//                            long pos = encoded.align(m);
                            long pos2[] = encoded.align2ComplimentV2(m);
//                            long pos2[] = encoded.fullAlign(m);
                            long pos = -1;
                            if(pos2!=null&&pos2.length>0){
                                pos = pos2[0];
                                pos = pos2.length;
                                //merMap = res.addResult(m, chr.getChrNumber(), pos2);
                            }
                            
                            int idx = (int) (pos-i);
                            if(pos<0){
                              idx = 0;
                            }
//                            System.out.println(""+chr.getName()+" "+sub+" "+sub.length()+" : "+m+" pos : "+pos+" idx : "+idx);
                            
                            System.out.print("\t"+pos);
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("Before Mermap this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            merMap = res.addResult(m, chr.getChrNumber(), pos2); // Still confuse After pss this line all member in pos2 will chang from 28 bit of position to 36 bit of chr:Pos Hoe did it happen !!!
                            // But it work!
//                            if(pos2!=null){
//                                    System.out.println();
//                                    System.out.println("After mer map this is pos2 Check: before add to mer read = " + pos2[0]);
//                                    System.out.println();
//                                }
//                            res.addResultV2(m, chr.getChrNumber(), pos2, seq.getReadName());
                            //System.out.println("Check seq length" + seq.getShortgunLength());
                            int totalMer = (seq.getShortgunLength()-mer)+1;
                            
                            /*************************************************************************************************************/
                            /* -------------------------New Implement Part (Stroe in object)---------------------------------------------*/
                            if (seq.getMerReadSize() < totalMer){
                                MerRead merRead = new MerRead();
                                
                                merRead.addMatchResultStrand(m, pos2, index, chr.getChrNumber());
                                seq.addMerRead(merRead);
//                                System.out.println("(First Time) Size Mer Read check: " + seq.getMerReadSize());
                            }else{
//                                System.out.println("(Other time) Size Mer Read check: " + seq.getMerReadSize());
                                if (index<seq.getMerReadSize()){
//                                    System.out.println("Check Index: " + index);
                                    MerRead merRead = seq.getMerRead().get(index); // call back existing merRead to update
                                    merRead.addMatchResultStrand(m, pos2, index, chr.getChrNumber());
                                    seq.addMerReadByIndex(index,merRead);
//                                    System.out.println("(Other time) Should be constant: Size Mer Read check: " + seq.getMerReadSize());
                                }                                   
                            }
                            /*-----------------------------------------------------------------------------------------------------------*/
                            /*************************************************************************************************************/
                            
                        }
                    //System.out.println(" This mer Map check: "+ (merMap == null));
                    res.createMap(seq.getReadName(), merMap);                    
                    }
                     System.out.println();
                    /* New Implement Part */
                    //seq.countAlignmentData(); // Create Alignment count data before change ShortgunSequence
                    
                    /*--------------------*/
                }
                /* End */
                
                encoded.lazyLoad();
                encoded = null;
                
                System.gc();
                
                
                
            } catch (IOException ex) {
                Logger.getLogger(BinaryAligner.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        //AlignmentResultRead alinResult = new AlignmentResultRead();  
        Enumeration<ShortgunSequence> seqs = input.getInputSequence().elements();
        while(seqs.hasMoreElements()){
                    ShortgunSequence seq = seqs.nextElement();
                    seq.countAlignmentData();
                    alinResult.addResult(seq);  
        }
        
        //return res;
        return alinResult;
    }
    
}
