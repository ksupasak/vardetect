/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import java.io.IOException;
import java.io.RandomAccessFile;

/**
 *
 * @author worawich
 */
public class TestRandomAccess {
    
    public static void main(String[] args){
        String fileName = "";
        
        try {
            // create a new RandomAccessFile with filename test
            RandomAccessFile raf = new RandomAccessFile(fileName, "rw");

            // write something in the file
            raf.writeUTF("Hello World");

            // set the file pointer at 0 position
            raf.seek(0);

            // print the string
            System.out.println("" + raf.readUTF());

            // set the file pointer at 5 position
            raf.seek(5);

            // write something in the file
            raf.writeUTF("This is an example");

            // set the file pointer at 0 position
            raf.seek(0);

            // print the string
            System.out.println("" + raf.readUTF());
         
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        
    }
    
}
