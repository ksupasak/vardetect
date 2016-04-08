/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 *
 * @author soup
 */
public class NGSCMD2 {
    public static void main(String args[]) throws FileNotFoundException, IOException{
        long i = 1231223123;
        DataOutputStream os = new DataOutputStream(new FileOutputStream("binout.dat"));
        os.writeLong(i);
        i++;
        os.writeLong(i);
        
        os.close();
        
        DataInputStream is = new DataInputStream(new FileInputStream("binout.dat"));
        long j = is.readLong();
        j = is.readLong();
        is.close();
        System.out.println(j);
        
        
    }
}
