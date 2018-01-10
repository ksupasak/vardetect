/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.cmd;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class FilterReadNameFromTargetGrep {
    public static void main(String[] args) throws FileNotFoundException, IOException{
        /**
         * Code below use to extract read name from variant report (this file is a file that we use grep command to cut out only interest group along with the first pattern)
         * and write only read name to new file
         */
        String file = "/Volumes/PromisePegasus/worawich/Download_dataset/Manon_TB/Martin_new_list/Martin_more/RD239.txt";
        ArrayList<String> readList = new ArrayList();
        RandomAccessFile rb = new RandomAccessFile(file,"rw");
        String line = "";
        while ((line = rb.readLine()) != null) {
            
            if(line.charAt(0)=='G'){
                
            }else if(line.charAt(0) == '-'){
            }else{
                String[] compLine = line.split(",");
                String[] readNameFull = compLine[9].split("\\.");
                String read = readNameFull[0];
                
                if(!readList.contains(read)){
                    readList.add(read);
                }
            }
        }
        
        FileWriter writer;
        /**
         * Check File existing
         */
        
        File f = new File(file.split("\\.")[0]+"_readName.txt"); //File object
        writer = new FileWriter(f);
        for(int i = 0;i<readList.size();i++){
            writer.write(readList.get(i));
            writer.write("\n");
        }
        
        writer.flush();
        writer.close();
        /*******************************************************/
    }
}
