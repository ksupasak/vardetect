/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import java.util.Random;

/**
 *
 * @author soup
 */
public class SimulatorUtil {
    
    
    
    
    public static InputSequence simulateIndel(ChromosomeSequence chr, int number, int shortgun){
        
        InputSequence is = new InputSequence();
        StringBuffer sb = chr.getSequence();
        Random rand = new Random();
        for(int i = 0;i<number;i++){
            
            int pos = rand.nextInt(sb.length()-2000)+1000;
            int type = rand.nextInt(2);
            int length = rand.nextInt(5)+1;
            
            String template;
            StringBuffer alt = new StringBuffer();
            if(type==0){ //insertion
                for(int j =0 ;j<length;j++)alt.append("A");
                template = sb.substring(pos-shortgun,pos)+alt+sb.substring(pos, pos+shortgun);
         
            }else{ //deletion
                template = sb.substring(pos-shortgun,pos)+sb.substring(pos+length, pos+shortgun);
            }
            
            if(template.indexOf("N")>=0)
                i--;
            else{
                System.out.println("Pos : " + pos+"\t"+"Type : " + type+"\t"+"Length : " + length);
                System.out.println(template);

                ShortgunSequence ss = new ShortgunSequence(template);
                is.addRead(ss);    
            } 
            
        
        }
        
        
        
        return is;
        
        
    }
}
