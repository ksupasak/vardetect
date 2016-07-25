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
public class SimulatorUtil_WholeGene {
    
    public static InputSequence simulateWholeGene(ReferenceSequence ref, int num_read, int ln_read, String numberchrA, String numberchrB){
        if (numberchrA == "24"){
            numberchrA = "Y";
        }else if(numberchrA == "23"){
            numberchrA = "X";
        }
        
        if (numberchrB == "24"){
            numberchrB = "Y";
        }else if(numberchrB == "23"){
            numberchrB = "X";
        }
        
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
        
            System.out.println("chr number: " + chrNum);
            ChromosomeSequence chr = chrs.elementAt(chrNum);
            System.out.println("Chromosome name: "+chr.getName());

            if (chr.getName().equalsIgnoreCase(namechrA)){
                chrA = chr;
            }
            else if (chr.getName().equalsIgnoreCase(namechrB)){
                chrB = chr;
            }   
        }
        System.out.println("concatenate");
        CharSequence iniTemplate = SequenceUtil.concatenateChromosome(chrA, chrB, ln_read-1, ln_read-1);
        
        String read;
        String readName;
        for(int i = 0;i<num_read;i++){
            int iniread =  0;
            int overLimitCheck = 0;
            while(overLimitCheck == 0){
                iniread = rand.nextInt(ln_read);
                if (iniread<(ln_read-1)){
                    overLimitCheck = 1;
                }
            }
            
            System.out.println("iniread : "+iniread);
            System.out.println("ln_read : "+ln_read);
            read = iniTemplate.subSequence(iniread, iniread+ln_read).toString();
            
            
            System.out.println("Initial position: " + iniread);
            System.out.println("Number of base from " + namechrA + " : " + ((ln_read-iniread)-1));
            System.out.println("Number of base from " + namechrB + " : " + (iniread+1));
            System.out.println(read);
            
            ShortgunSequence ss = new ShortgunSequence(read);
            readName = "Read"+i;
            ss.addReadName(readName);
            is.addRead(ss);
            
        }
        
        //is.addName(ref.getChromosomes());
        
        
        //while(chrs.elements()..hasMoreElements()){
            
            //ChromosomeSequence chr = chrs.elements().nextElement();
            
            //System.out.println(chr.getName());
            
        //}
        /*System.out.println(chrs.elements().hasMoreElements());*/
        
        
        
        
        
        
        
        return is;
    }
    
    public static InputSequence simulateWholeGeneRandom(ReferenceSequence ref, int num_read, int ln_read, int num_shortgun){
        /* num_read is number of read per random time */
        InputSequence is = new InputSequence();
        Random rand1 = new Random(); /* For random chromosome */
        for(int i = 0;i<num_read;i++){
            String numberChrA = Integer.toString(rand1.nextInt(25-1)+1);
            String numberChrB = Integer.toString(rand1.nextInt(25-1)+1);

            //InputSequence is = simulateWholeGene(ref,num_shortgun,ln_read,Integer.toString(numberChrA),Integer.toString(numberChrB));
            
       
        
            System.out.println("Begin Simulate");
            
            if (numberChrA.equalsIgnoreCase("24")){
                numberChrA = "Y";
            }else if(numberChrA.equalsIgnoreCase("23")){
                numberChrA = "X";
            }

            if (numberChrB.equalsIgnoreCase("24")){
                numberChrB = "Y";
            }else if(numberChrB.equalsIgnoreCase("23")){
                numberChrB = "X";
            }
            
            String namechrA = "chr"+ numberChrA;
            String namechrB = "chr"+ numberChrB;
            ChromosomeSequence chrA = null,chrB = null;
            
            Random rand2 = new Random(); /* For random positionn on cancatenate sequence */
        
            Vector<ChromosomeSequence> chrs = ref.getChromosomes();
            System.out.println(chrs.size());
            System.out.println("Chromosome loop");
            System.out.println("namechrA => " + namechrA);
            System.out.println("namechrB => " + namechrB);
            for(int chrNum=0;chrNum<chrs.size();chrNum++){

                System.out.println("chr number: " + chrNum);
                ChromosomeSequence chr = chrs.elementAt(chrNum);
                System.out.println("Chromosome name: "+chr.getName());

                if (chr.getName().equalsIgnoreCase(namechrA) && chr.getName().equalsIgnoreCase(namechrB)){
                    chrA = chr;
                    chrB = chr;
                }else{
                    if(chr.getName().equalsIgnoreCase(namechrA)){
                        chrA = chr;
                    }else if(chr.getName().equalsIgnoreCase(namechrB)){
                        chrB = chr;
                    } 
                }
            }
            System.out.println("concatenate");
            CharSequence iniTemplate = SequenceUtil.concatenateChromosome(chrA, chrB, ln_read-1, ln_read-1);
            
            chrA.lazyLoad();
            chrB.lazyLoad();
            
            String read;
            String readName;
            for(int j = 0;j<num_shortgun;j++){
                int iniread =  0;
                int overLimitCheck = 0;
                while(overLimitCheck == 0){
                    iniread = rand2.nextInt(ln_read);
                    if (iniread<(ln_read-1)){
                        overLimitCheck = 1;
                    }
                }

                System.out.println("iniread : "+iniread);
                System.out.println("ln_read : "+ln_read);
                read = iniTemplate.subSequence(iniread, iniread+ln_read).toString();


                System.out.println("Initial position: " + iniread);
                System.out.println("Number of base from " + namechrA + " : " + ((ln_read-iniread)-1));
                System.out.println("Number of base from " + namechrB + " : " + (iniread+1));
                System.out.println(read);

                ShortgunSequence ss = new ShortgunSequence(read);
                readName = "Read"+i+"SS"+j;
                ss.addReadName(readName);
                is.addRead(ss);
                ss = null;
                chrs =null;
                
                System.gc();
            }
           
        }
        
        //is.addName(ref.getChromosomes());
        
        
        //while(chrs.elements()..hasMoreElements()){
            
            //ChromosomeSequence chr = chrs.elements().nextElement();
            
            //System.out.println(chr.getName());
            
        //}
        /*System.out.println(chrs.elements().hasMoreElements());*/
        
        
        
        
        
               
        return is;
    }
    
    public static InputSequence simulateComplexWholeGeneRandom(ReferenceSequence ref, int num_read, int ln_read, int num_shortgun){
        /* num_read is number of read per random time */
        InputSequence is = new InputSequence();
        Random rand1 = new Random(); /* For random chromosome */
        for(int i = 0;i<num_read;i++){
            String numberChrA = Integer.toString(rand1.nextInt(25-1)+1);
            String numberChrB = Integer.toString(rand1.nextInt(25-1)+1);

            //InputSequence is = simulateWholeGene(ref,num_shortgun,ln_read,Integer.toString(numberChrA),Integer.toString(numberChrB));
            
       
        
            System.out.println("Begin Simulate");
            
            if (numberChrA.equalsIgnoreCase("24")){
                numberChrA = "Y";
            }else if(numberChrA.equalsIgnoreCase("23")){
                numberChrA = "X";
            }

            if (numberChrB.equalsIgnoreCase("24")){
                numberChrB = "Y";
            }else if(numberChrB.equalsIgnoreCase("23")){
                numberChrB = "X";
            }
            
            String namechrA = "chr"+ numberChrA;
            String namechrB = "chr"+ numberChrB;
            ChromosomeSequence chrA = null,chrB = null;
            
            Random rand2 = new Random(); /* For random positionn on cancatenate sequence */
        
            Vector<ChromosomeSequence> chrs = ref.getChromosomes();
            System.out.println(chrs.size());
            System.out.println("Chromosome loop");
            System.out.println("namechrA => " + namechrA);
            System.out.println("namechrB => " + namechrB);
            for(int chrNum=0;chrNum<chrs.size();chrNum++){

                System.out.println("chr number: " + chrNum);
                ChromosomeSequence chr = chrs.elementAt(chrNum);
                System.out.println("Chromosome name: "+chr.getName());

                if (chr.getName().equalsIgnoreCase(namechrA) && chr.getName().equalsIgnoreCase(namechrB)){
                    chrA = chr;
                    chrB = chr;
                }else{
                    if(chr.getName().equalsIgnoreCase(namechrA)){
                        chrA = chr;
                    }else if(chr.getName().equalsIgnoreCase(namechrB)){
                        chrB = chr;
                    } 
                }
            }
            System.out.println("concatenate");
            CharSequence iniTemplate = SequenceUtil.concatenateComplexChromosome(chrA, chrB, ln_read-1, ln_read-1);
            
            chrA.lazyLoad();
            chrB.lazyLoad();
            
            String read;
            String readName;
            /* Loop for shortgun */
            for(int j = 0;j<num_shortgun;j++){
                int iniread =  0;
                int overLimitCheck = 0;
                while(overLimitCheck == 0){
                    iniread = rand2.nextInt(ln_read);
                    if (iniread<(ln_read-1)){
                        overLimitCheck = 1;
                    }
                }

                System.out.println("iniread : "+iniread);
                System.out.println("ln_read : "+ln_read);
                read = iniTemplate.subSequence(iniread, iniread+ln_read).toString();


                System.out.println("Initial position: " + iniread);
                System.out.println("Number of base from " + namechrA + " : " + ((ln_read-iniread)-1));
                System.out.println("Number of base from " + namechrB + " : " + (iniread+1));
                System.out.println(read);

                ShortgunSequence ss = new ShortgunSequence(read);
                readName = "Read"+i+"SS"+j;
                ss.addReadName(readName);
                is.addRead(ss);
                ss = null;
                chrs =null;
                
                System.gc();
            }
           
        }
        
        //is.addName(ref.getChromosomes());
        
        
        //while(chrs.elements()..hasMoreElements()){
            
            //ChromosomeSequence chr = chrs.elements().nextElement();
            
            //System.out.println(chr.getName());
            
        //}
        /*System.out.println(chrs.elements().hasMoreElements());*/
        
        
        
        
        
               
        return is;
    }
}
