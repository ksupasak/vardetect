/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
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
    Map<String,Integer> chrIndex = new LinkedHashMap();                 // store chr Index (Map that link btw original chr name and new chr name)
    ArrayList<Long> annoBinaryTree = new ArrayList();                   // array of long of startCode and stopCode (sorted) [chr 5bit][start position 28bit][annotation index 31bit]
    Map<Integer,Annotation> annoIndex = new LinkedHashMap();            // store Anootation and it index
 
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
         * this TreeMap for matching purpose
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
    
    public void putChrIndex(Map<String,Integer> inChrIndex){
        this.chrIndex = inChrIndex;
    }
    
    public void putAnnotationBinaryTree(ArrayList<Long> inAnnoBinaryTree){
        this.annoBinaryTree = inAnnoBinaryTree;
    }
    
    public void putAnnotationIndex(Map<Integer,Annotation> inAnnoIndex){
        this.annoIndex = inAnnoIndex;
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
