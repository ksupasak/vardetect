/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.AlignmentResult;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author worawich
 */
public class VisualizeResult {
    
    
    public static void visualizeAlignmentResult(AlignmentResult inRes){
        Map<String,ArrayList<Map>> frontMapRes = inRes.getAlignmentResult();
        
        Set set = frontMapRes.keySet();
        Iterator iter = set.iterator();
        
        while (iter.hasNext()){
            Object rdName = iter.next();
            ArrayList<Map> merMap = frontMapRes.get(rdName);                        
            System.out.println("\nThis is result check(read name) : "+rdName);
                        
            long array_size = merMap.size();
            //System.out.println("This is size check for number of mer in current read size is : " + array_size);
            Map<Long,long[]> subMap = new HashMap();
            for(int i = 0;i<array_size;i++){
                subMap = merMap.get(i);
                //System.out.println("Check subMap is emty? : "+subMap.isEmpty());
                Set setMer = subMap.keySet();
                Iterator iterMer = setMer.iterator();
                //System.out.println("check Key size: "+setMer.size());
                if(subMap.isEmpty()){
                    
                    //System.out.println("subMap is empty");
                }else{
                    while(iterMer.hasNext()){
                        System.out.print("\n");
                        Object dum = iterMer.next();
                        //System.out.println("Key is: "+dum);
                        long[] codePos = subMap.get(dum);
                        //System.out.println("CodePos: "+ codePos + " "+codePos.length);
                        //System.out.println("YOYO CheckCheck");
                        System.out.print("Mer code: "+ dum);

                        for (int j=0;j<codePos.length;j++){

                            long chrnumber = codePos[j]>>28;
                            long position = codePos[j]&268435455;

                            System.out.println("\tchr: "+chrnumber + "\tPosition: "+position);
                        }
                    
                    
                    }
                }
                
            }
            
            
        }
        //SefrontMapRes
        
        //frontMapRes.get(inRes);
    }
    
}
