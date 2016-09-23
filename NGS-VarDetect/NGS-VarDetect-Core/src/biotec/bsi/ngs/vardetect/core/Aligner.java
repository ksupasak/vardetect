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
    public AlignmentResultRead align(ReferenceSequence ref, InputSequence input);
    public AlignmentResultRead align(InputSequence input);
    public AlignmentResultRead alignV2(ReferenceSequence ref, InputSequence input);
    public AlignmentResultRead alignV2(InputSequence input);
    public AlignmentResultRead alignMultithread(ReferenceSequence ref, InputSequence input, int numThread)throws InterruptedException;
    public AlignmentResultRead alignMultithread(InputSequence input, int numThread)throws InterruptedException;
    
}
