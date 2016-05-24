/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.alignment;

import biotec.bsi.ngs.vardetect.core.Aligner;

/**
 *
 * @author soup
 */
public class AlignerFactory {
    public static Aligner getAligner(){
        return new BinaryAligner();
    }
}
