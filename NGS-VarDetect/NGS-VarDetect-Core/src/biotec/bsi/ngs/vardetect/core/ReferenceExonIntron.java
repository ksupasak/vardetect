/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.Collections;
import java.util.Vector;

/**
 *
 * @author Worawich
 */
public class ReferenceExonIntron {
    
    Vector<ExonIntron> data;
    
    public ReferenceExonIntron(){
        data = new Vector<ExonIntron>();
    }
    
    public void addData(ExonIntron in){
        data.add(in);
    }
    public void reverseorder(){
        Collections.reverse(data);
    }
     public Vector<ExonIntron> getdata(){
        return data;
    }
}
