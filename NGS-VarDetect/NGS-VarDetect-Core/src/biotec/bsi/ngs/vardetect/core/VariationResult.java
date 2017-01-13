/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author worawich
 */
public class VariationResult {
    private Map<Integer,String[]> variationResult;
    
    public VariationResult(){
        this.variationResult = new LinkedHashMap();
        
    }
    
    public void addVariationMap(Map<Integer,String[]> variation){
       this.variationResult.putAll(variation);
    }
            
    
}
