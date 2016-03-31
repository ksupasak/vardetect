/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;
import java.io.FileNotFoundException;
import java.io.PrintStream;
/**
 *
 * @author soup
 */
public class ChromosomeSequence {
    String name;
    StringBuffer seq;
    public ChromosomeSequence(String name, StringBuffer seq){
        this.name = name;
        this.seq = seq;
    }

    public void writeToPath(String path, String fa) throws FileNotFoundException {

       PrintStream ps = new PrintStream(path+".fa");
       ps.println(">"+name);
       ps.println(seq);

    }
    
}
