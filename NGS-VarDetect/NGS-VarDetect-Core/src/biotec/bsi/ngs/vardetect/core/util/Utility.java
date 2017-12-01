/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 *
 * @author worawich
 */
public class Utility {
    
    public static void modifyFile(String file) throws IOException{
        /**
         * read file and modify then write to new file
         */
        
        FileWriter writer;
        String[] dummysaveFile = file.split("\\.");
        String saveFile = dummysaveFile[0] + "_Modif." + dummysaveFile[1];
        /**
         * Check File existing
         */
        
        File f = new File(saveFile); //File object        
        if(f.exists()){
            // append if exist
//            ps = new PrintStream(new FileOutputStream(filename,true));
            writer = new FileWriter(saveFile,true);
        }else{
            // create new
//            ps = new PrintStream(filename);
            writer = new FileWriter(saveFile);
        }
        writer = new FileWriter(saveFile);
        
        Charset charset = Charset.forName("US-ASCII");
        Path path = Paths.get(file);
        
        StringBuffer seq = new StringBuffer();

        try (BufferedReader reader = Files.newBufferedReader(path, charset)) {
            String line = null;
            String[] data = null; 
            while ((line = reader.readLine()) != null) {
                data = line.split("_");
                writer.write(data[0]);
                writer.write("\n");
            }
        }
        
        writer.flush();
        writer.close();
    } 
}
