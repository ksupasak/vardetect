/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import java.io.EOFException;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 *
 * @author worawich
 */
public class TestRandomAccess {
    
    public static void main(String[] args){
        String fileName = "/Volumes/PromisePegasus/worawich/Referense/TB_reference/H37Rv_NC_000962_reIndex.fa";
        
        try {
            // create a new RandomAccessFile with filename test
            RandomAccessFile raf = new RandomAccessFile(fileName, "rw");
            // set the file pointer at 0 position
            raf.seek(10149);
            
            // print the string
            System.out.println("Test random access");
            System.out.println("" + raf.readLine());
            System.out.println("Pointer is on : " + raf.getFilePointer());
            System.out.println("" + raf.readLine());
            System.out.println("Pointer is on : " + raf.getFilePointer());
            System.out.println("" + raf.readLine());
       
        } catch (EOFException e) {
            
        } catch (IOException ex) {
            ex.printStackTrace();
        }
   
    }
    
}
