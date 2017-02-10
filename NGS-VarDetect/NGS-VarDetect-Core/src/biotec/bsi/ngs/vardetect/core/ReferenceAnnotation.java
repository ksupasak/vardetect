/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

/**
 *
 * @author Worawich
 */
public class ReferenceAnnotation {
    
    ArrayList<Annotation> data;
    Map<Long,Map<Long,ArrayList<Annotation>>> refAnno; 
    
    public ReferenceAnnotation(){
        data = new ArrayList<Annotation>();
        refAnno = new TreeMap();
    }
    
    public void putData(Annotation inAnno,long startCode,long stopCode){
        /**
         * put annotation object data into TreeMap
         * TreeMap has two layer 
         *  1. key is startCode : value is TreeMap
         *  2. key is stopCode : value is ArrayList<Annotation>
         * this TreeMap for matching propose
         */
        if(this.refAnno.containsKey(startCode)){
            Map<Long,ArrayList<Annotation>> refII = refAnno.get(startCode);
            if(refII.containsKey(stopCode)){
                ArrayList<Annotation> listAnno = refII.get(stopCode);
                listAnno.add(inAnno);
                refII.put(stopCode, listAnno);
            }else{
                ArrayList<Annotation> listAnno = new ArrayList();
                listAnno.add(inAnno);
                refII.put(stopCode, listAnno);
            }
            this.refAnno.put(startCode, refII);
        }else{
            Map<Long,ArrayList<Annotation>> refII = new TreeMap();
            ArrayList<Annotation> listAnno = new ArrayList();
            listAnno.add(inAnno);
            refII.put(stopCode, listAnno);
            this.refAnno.put(startCode, refII);
        }
    }
    
    public void addData(Annotation in){
        data.add(in);
    }
    public void reverseorder(){
        Collections.reverse(data);
    }
    public ArrayList<Annotation> getdata(){
        return data;
    }
    public Map<Long,Map<Long,ArrayList<Annotation>>> getReferenceAnnotation(){
        return this.refAnno;
    }
}
