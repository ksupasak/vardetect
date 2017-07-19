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
    Map<String,Long> chrIndex = new LinkedHashMap();                 // store chr Index (Map that link btw original chr name and new chr name)
    Map<Long,String> chrIndexReverse = new LinkedHashMap();
    ArrayList<Long> annoBinaryTree = new ArrayList();                   // array of long of startCode and stopCode (sorted) [chr 5bit][start position 28bit][annotation index 31bit]
    Map<Integer,Annotation> annoIndex = new LinkedHashMap();            // store Anootation and it index
    private long indexAnnoMask = 8388607;                               // mask of 23bit
    private long chrPosCompareMask = -8388608;                          // use to tranform chrPosAnnoIndex code into chrPos (bit 1) and bit 0 at AnnoIndex Ex. 111100000 (front part have value but back part is not) for easy binary search comparing
    
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
    
    public void putChrIndex(Map<String,Long> inChrIndex,Map<Long,String> inChrIndexReverse){
        this.chrIndex = inChrIndex;
        this.chrIndexReverse = inChrIndexReverse;
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
    public Map<Integer,Annotation> getAnnotationIndex(){
        return this.annoIndex;
    }
    public Map<String,Long> getChrIndex(){
        return this.chrIndex;
    }
    public Map<Long,String> getChrIndexReverse(){
        return this.chrIndexReverse;
    }
    
    public int mapToAnotationBinaryTreeWithPosStart(long chrPosStart, long chrPosStop){
        /**
         * This function try to map chrPosStart With Binary search to the Array of chrPosAnnoIndex
         * then return the Annotation index that this chrPosStart has match or in the rage of start and stop point of that annotation
         * 
         * If chrPosStart did not match to any annotation range the program will move to use chrPosStop instead
         * If chrPosStop not match, we will return -1 (this match pattern is not in the range of any annotation)
         * 
         */
        
        int iniPoint = 0;
        int lastPoint = annoBinaryTree.size();
        int midPoint = (lastPoint + iniPoint)/2;
        long compareChrPos;
        int annotationIndex = 0;
        char status = 'n';              // status indicate that the input is lower or higher than compareChrPos 'l' is lower |  'h' is higher (if 'l' it mean the chrPosStart is lower than that mid point)
        
        while(true){
            compareChrPos = annoBinaryTree.get(midPoint) & this.chrPosCompareMask;
            
            if(compareChrPos == chrPosStart){
                annotationIndex = (int)(this.annoBinaryTree.get(midPoint)&this.indexAnnoMask);
                return annotationIndex;
            }else if(chrPosStart < compareChrPos){
                status = 'l';
                lastPoint = midPoint-1;
                if(iniPoint>lastPoint){
                    break;
                }else{
                    midPoint = (lastPoint + iniPoint)/2;
                }
            }else if(chrPosStart > compareChrPos){
                status = 'h';
                iniPoint = midPoint+1;
                if(iniPoint>lastPoint){
                    break;
                }else{
                    midPoint = (lastPoint + iniPoint)/2;
                }               
            }
            
//            if(iniPoint>lastPoint){
//                // has case check to compensate to the real mid point that the chrPos has lower or higher. because midpoint has being update to new midpoint and then the while loop is break before it has been use to compare
//                // So, the status is stiil idicate lower or higher which does not mean this new midpoint but mean the midpoint in the past
//                if(status == 'l'){
//                    midPoint++;
//                }else if(status == 'h'){
//                    midPoint--;
//                }
//                break;
//            }
        }
        
        if(status=='l'){
            if(midPoint == 0){
                
                /**
                 * midPoint is Start Point
                 * it's lower than start
                 * it have a chance that back part of it may overlap with this annotation. So, we have to check with chrPosStop.
                 * So we compare chrPosstop with chrPos of this annotation (midPoint)
                 */
                compareChrPos = annoBinaryTree.get(midPoint) & this.chrPosCompareMask;
                if(chrPosStop < compareChrPos){
                    // chrPosStop has lower than compareChrPos. This mean back part is not in the range of this annotation
                    return -1;
                }else if(chrPosStop>=compareChrPos){
                    // chrPosStop is higher than compareChrPos this mean it's back part is overlap in the range of this annotation. So we return this anootation index.
                    int midPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint) & this.indexAnnoMask);
                    annotationIndex = midPointAnnotationIndex;
                    return annotationIndex;
                }
                
                return -1;

            }else if(midPoint == annoBinaryTree.size()){
                int midPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint) & this.indexAnnoMask);
                /**
                * midpoint is Stop Point
                * It's lower than stop
                */
                annotationIndex = midPointAnnotationIndex;
                return annotationIndex;
                
            }else{
            
                int midPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint) & this.indexAnnoMask);
                int leftPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint-1) & this.indexAnnoMask);
                int rightPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint+1) & this.indexAnnoMask);

                if(leftPointAnnotationIndex == midPointAnnotationIndex){
                    /**
                     * midpoint is Stop Point
                     * It lower than stop
                     */

                    annotationIndex = midPointAnnotationIndex;
                    return annotationIndex;

                }else if(rightPointAnnotationIndex == midPointAnnotationIndex){
                    /**
                     * midPoint is Start Point
                     * it's lower than start
                     * it have a chance that back part of it may overlap with this annotation. So, we have to check with chrPosStop.
                     * So we compare chrPosstop with chrPos of this annotation (midPoint)
                     */
                    compareChrPos = annoBinaryTree.get(midPoint) & this.chrPosCompareMask;
                    if(chrPosStop < compareChrPos){
                        // chrPosStop has lower than compareChrPos. This mean back part is not in the range of this annotation
                        return -1;
                    }else if(chrPosStop>=compareChrPos){
                        // chrPosStop is higher than compareChrPos this mean it's back part is overlap in the range of this annotation. So we return this anootation index.
                        annotationIndex = midPointAnnotationIndex;
                        return annotationIndex;
                    }

                }
            }
        }else if(status=='h'){
            
            if(midPoint == 0){
                int midPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint) & this.indexAnnoMask);
                /**
                * midPoint is Start point
                * it's higher than start
                */
                
                annotationIndex = midPointAnnotationIndex;
                return annotationIndex;

            }else if(midPoint == annoBinaryTree.size()){
                /**
                * midpoint is Stop Point
                * It's higher than stop
                */
                
                return -1;
                
            }else{

                int midPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint) & this.indexAnnoMask);
                int leftPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint-1) & this.indexAnnoMask);
                int rightPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint+1) & this.indexAnnoMask);

                if(leftPointAnnotationIndex == midPointAnnotationIndex){
                    /**
                     * midPoint is Stop point
                     * It's higher than stop
                     * it have a chance that back part of it may overlap with next annotation. So, we have to check with chrPosStop.
                     * So we compare chrPosstop with chrPos of next annotation (next anno start point)
                     */
                    
                    compareChrPos = annoBinaryTree.get(midPoint+1) & this.chrPosCompareMask;
                    if(chrPosStop < compareChrPos){
                        // chrPosStop has lower than or equal to compareChrPos. This mean back part is not in the range of next annotation index
                        return -1;
                    }else if(chrPosStop>=compareChrPos){
                        // chrPosStop higher or equal with compareChrPos. this mean back part is overlap with next annotation
                        midPointAnnotationIndex = (int)(this.annoBinaryTree.get(midPoint+1) & this.indexAnnoMask);
                        annotationIndex = midPointAnnotationIndex;
                        return annotationIndex;
                    }
                    
                    annotationIndex = midPointAnnotationIndex;
                    return annotationIndex;
                }else if(rightPointAnnotationIndex == midPointAnnotationIndex){
                    /**
                     * midPoint is Start point
                     * it's higher than start
                     */
                    
                    annotationIndex = midPointAnnotationIndex;
                    return annotationIndex;
                }
            }
        }
        
//        check correction of this logic 
        return -1;
    }
    
    public void mapToAnotationBinaryTreeWithPosStop(long chrPosStop){
        
    }
}
