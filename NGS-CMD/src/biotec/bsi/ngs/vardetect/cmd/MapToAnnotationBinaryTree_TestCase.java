/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.ReferenceAnnotation;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author worawich
 */
public class MapToAnnotationBinaryTree_TestCase {
    public static void main(String[] args) throws IOException {
        String answer = "";
        String key = "11236544"; 
        /**
         * Generate reference annotation
         */
        ArrayList<Long> inputList = new ArrayList();
        ArrayList<Long> annotation = new ArrayList();
        long[] chrList = {1,2};
        long[] startstop = {10,20,500,1000,1005,1010};
        int count = 0;
        int annoIndex = 1;
         
        for(int i =0;i<chrList.length;i++){
            for(int j =0;j<startstop.length;j++){
                count++;
                long chr = chrList[i];
                long pos = startstop[j];
                
                long chrPosAnnoIndex = (((chr<<28)+pos)<<23)+annoIndex;
                annotation.add(chrPosAnnoIndex);
                
                if(count==2){
                    count=0;
                    annoIndex++;
                } 
            }
        }
        Collections.sort(annotation);
        
        ReferenceAnnotation refAnno = new ReferenceAnnotation();
        refAnno.putAnnotationBinaryTree(annotation);
        
        /**
         * Generate input
         */
        long[] inputChrList = {1,2};
        long[] inputStartStop = {5,10,0,15,700,1001,1009,1500,1001,1006,30,600,10,15,5,10};
        count=0;
        int chrindex =0;
        for(int j =0;j<inputStartStop.length;j++){
            count++;
            long chr = inputChrList[chrindex];
            long pos = inputStartStop[j];

            long chrPosAnnoIndex = ((chr<<28)+pos)<<23;
            inputList.add(chrPosAnnoIndex);

            if(count%8==0){
                chrindex++;
            } 
        }
        
        /**
         * test map to annotation reference 
         */
        
        
        for(int i = 0;i<inputList.size();i=i+2){
            long chrPosStart = inputList.get(i);
            long chrPosStop = inputList.get(i+1);
            
            int returnIndex = refAnno.mapToAnotationBinaryTreeWithPosStart(chrPosStart, chrPosStop);
            answer = answer + String.valueOf(returnIndex);
        }
        
        
        if(answer.endsWith(key)){
            System.out.println("Correct");
        }else{
            System.out.println("no correct");
            System.out.println("the result is :" + answer);
        }
        System.out.println();
    }
}
