/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.ChromosomeSequence;
import biotec.bsi.ngs.vardetect.core.EncodedSequence;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.MapResult;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.ShortgunSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_aon.*;
import static biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_aon.simulateWholeGene;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Worawich
 */
public class NGSCMD3 {
    
    public static void main(String[] args) throws IOException{
        
        String file_path = "/Users/worawich/VMdev/dataScieneToolBox/projects/genomic_projects/Reference_hg19/Map_Result.map";
        
        
        
        //ReferenceSequence ref = SequenceUtil.readReferenceSequence(args[1]);
        ReferenceSequence ref_read = SequenceUtil.readReferenceSequence(args[0]);
        System.out.println("Read Done");
        System.out.println("File path is : "+ ref_read.getPath());


        //////////ReferenceSequence ref = SequenceUtil.readReferenceSequence(args[1]);
        ///ChromosomeSequence chr = ref.getChromosomes().elementAt(0);
        ///System.out.println("Name of Pick chromosome : " + chr.getName());


        ///EncodedSequence encode = SequenceUtil.getEncodeSequence(chr);
        ///System.out.println("Size of encode referecce : " + encode.getEncodeMap().size());


        InputSequence is = simulateWholeGene(ref_read,5,100,20,21); 

        //Enumeration<ShortgunSequence> e = is.seqs.elements();
        MapResult mapRes = new MapResult();
        //Map<Long,Long> mapRes = new HashMap();
        //System.out.println("all key "+encode.getEncodeMap().keySet());
        
        mapRes = SequenceUtil.mapGenomeShotgunV3(ref_read, is,12,13); // Last 2 number is defined for element in reference that we want to map to. chr 20 and 21 is on element 12 ans 13
        /////mapRes.readFromPath(file_path, "map");
        
        //mapRes.writeToPath(ref.getPath(), "map");
        //System.out.println("Summary : All key contains = "+mapRes.getResultMap().keySet());
        System.out.println(mapRes.getResultArray().size());

        for (int i = 0;i<mapRes.getAlignPosition().size();i++){
            //dummy = mapRes.getResultArray().get(i);
            
            System.out.println("Read name : " + mapRes.getReadName().get(i) + " Align at : " + mapRes.getAlignPosition().get(i) + " : On chromosome : "+ mapRes.getchrNumber().get(i)+ " : Match : " + mapRes.getNumMatch().get(i));

        }

        SequenceUtil.createHistrogram(mapRes,"Read0");
        SequenceUtil.createHistrogram(mapRes,"Read1");
        SequenceUtil.createHistrogram(mapRes,"Read2");
        SequenceUtil.createHistrogram(mapRes,"Read3");
        SequenceUtil.createHistrogram(mapRes,"Read4");
        /*int count = 0;
        while(e.hasMoreElements()){
            System.out.println("Read Number : " + count++);
            ShortgunSequence ss = e.nextElement();

            EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(ss.seq);

            SequenceUtil.mapGenome(encode, encodeSim); // encode = referance map ; encodeSim = read

        }/*
        //Enumeration<ShortgunSequence> e = is.seqs.elements();


        /*
        ChromosomeSequence chr = ref.getChromosomes().elementAt(0);
        EncodedSequence encode = SequenceUtil.encodeSerialChromosomeSequence(chr);

        InputSequence is = SimulatorUtil.simulateIndel(chr, 5, 100);

        Enumeration<ShortgunSequence> e = is.seqs.elements();

        while(e.hasMoreElements()){

            ShortgunSequence ss = e.nextElement();

            EncodedSequence encodeSim = SequenceUtil.encodeSerialReadSequence(ss.seq);

            SequenceUtil.mapGenome(encode, encodeSim);


        }*/
        
    }
    
}
