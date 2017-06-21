/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

/**
 *
 * @author soup
 */
public interface Aligner {
    
    public void setReferenceSequence(ReferenceSequence ref);
    public AlignmentResultRead align(ReferenceSequence ref, InputSequence input, int numMer, int threshold);
    public AlignmentResultRead align(InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV2(ReferenceSequence ref, InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV2(InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV3(ReferenceSequence ref, InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV3(InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV4(ReferenceSequence ref, InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV4(InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV5(ReferenceSequence ref, InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignV5(InputSequence input, int numMer, int threshold);
    public AlignmentResultRead alignMultithread(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithread(InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadV3(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadV3(InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadV4(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadV4(InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadV5(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadV5(InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadLongRead(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadLongRead(InputSequence input, int numThread, int numMer, int threshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadLongReadRepeatCut(ReferenceSequence ref, InputSequence input, int numThread, int numMer, int threshold, int repeatThreshold)throws InterruptedException;
    public AlignmentResultRead alignMultithreadLongReadRepeatCut(InputSequence input, int numThread, int numMer, int threshold, int repeatThreshold)throws InterruptedException;
    public AlignmentResultRead localAlign(EncodedSequence ref, InputSequence input, int kmer, long numberOflocalRef);
    
}
