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
    public AlignmentResult align(ReferenceSequence ref, InputSequence input);
    public AlignmentResult align(InputSequence input);
    
    
}