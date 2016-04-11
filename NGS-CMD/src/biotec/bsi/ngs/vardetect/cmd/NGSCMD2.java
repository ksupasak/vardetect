/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import biotec.bsi.ngs.vardetect.core.ExonIntron;
import biotec.bsi.ngs.vardetect.core.InputSequence;
import biotec.bsi.ngs.vardetect.core.ReferenceExonIntron;
import biotec.bsi.ngs.vardetect.core.ReferenceSequence;
import biotec.bsi.ngs.vardetect.core.util.SequenceUtil;
import biotec.bsi.ngs.vardetect.core.util.SimulatorUtil_aon;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Random;
import java.util.Vector;

/**
 *
 * @author soup
 */
public class NGSCMD2 {
    
    public static void main(String args[]) throws FileNotFoundException, IOException{
        
        ReferenceExonIntron ee = SequenceUtil.readExonIntron(args[3]);
        Vector<ExonIntron> data = ee.getdata();
        System.out.println("Vector size" + data.size());
        System.out.print("Pick" + data.elementAt(500).getdirection());
        int iniguess = Math.abs(data.elementAt(500).getdirection());
        //for (int i == data.elementAt(500).getdirection())
        
        //ReferenceSequence ref = SequenceUtil.readReferenceSequence(args[1]);
        //InputSequence is = new InputSequence();  
        
        //is = SimulatorUtil_aon.simulateWholeGene(ref, 5, 100);
        
        //SequenceUtil.extractReferenceSequence(args[1]);
    }
}
