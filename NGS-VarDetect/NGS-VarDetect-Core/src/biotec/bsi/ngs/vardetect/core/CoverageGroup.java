/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class CoverageGroup {
    
    private ArrayList<Variation> listVariation;
    
    public CoverageGroup(Variation inVar){
        this.listVariation = new ArrayList();
        this.listVariation.add(inVar);  
    }
    
    public void addVariation(Variation inVar){
        this.listVariation.add(inVar);
    }
    
}
