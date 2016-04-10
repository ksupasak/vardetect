/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import java.util.Random;
import java.util.Vector;

/**
 *
 * @author Worawich
 */
public class SimulatorUtil_aon {
    
    public static InputSequence simulateWholeGene(ReferenceSequence ref, int num_read, int ln_read, int numberchrA,int numberchrB){
        System.out.println("Begin Simulate");
        String namechrA = "chr"+numberchrA;
        String namechrB = "chr"+numberchrB;
        ChromosomeSequence chrA = null,chrB = null;
        InputSequence is = new InputSequence();
        Random rand = new Random();
        
        Vector<ChromosomeSequence> chrs = ref.getChromosomes();
        System.out.println(chrs.size());
        System.out.println("Chromosome loop");
        for(int chrNum=0;chrNum<chrs.size();chrNum++){
        
            System.out.println("chr number" + chrNum);
            ChromosomeSequence chr = chrs.elementAt(chrNum);
            System.out.println("Chromosome name"+chr.getName());

            if (chr.getName().equalsIgnoreCase(namechrA)){
                chrA = chr;
            }
            else if (chr.getName().equalsIgnoreCase(namechrB)){
                chrB = chr;
            }   
        }
        System.out.println("boncatenate");
        CharSequence iniTemplate = SequenceUtil.concatenateChromosome(chrA, chrB, ln_read-1, ln_read-1);
        
        String read;
        for(int i = 0;i<num_read;i++){
            int iniread = rand.nextInt(ln_read);
            
            read = iniTemplate.subSequence(iniread, iniread+ln_read).toString();
            
            System.out.println("Initial position: " + iniread);
            System.out.println("Number of base from " + namechrA + " : " + ((ln_read-iniread)-1));
            System.out.println("Number of base from " + namechrB + " : " + (iniread+1));
            System.out.println(read);
            
            ShortgunSequence ss = new ShortgunSequence(read);
            is.seqs.add(ss);         
        }
        
        //while(chrs.elements()..hasMoreElements()){
            
            //ChromosomeSequence chr = chrs.elements().nextElement();
            
            //System.out.println(chr.getName());
            
        //}
        /*System.out.println(chrs.elements().hasMoreElements());*/
        
        
        
        
        
        
        
        return is;
    }
    
}
