/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
/**
 *
 * @author soup
 */
public class ChromosomeSequence {
    
    
    ReferenceSequence ref;
    
    String name;
    StringBuffer seq;
    Long chrNum;
    
    
    public ChromosomeSequence(ReferenceSequence ref, String name, StringBuffer seq){
        this.ref = ref;
        this.name = name;
        this.seq = seq;
    }

    public void writeToPath(String path, String fa) throws FileNotFoundException {

       PrintStream ps = new PrintStream(path+".fa");
       ps.println(">"+name);
       ps.println(seq);

    }
    
    
    public String getName(){
        return name;
    }
    
    public Long getChrNumber(){
        String[] dummy = name.split("r");
        if (dummy[1]=="X"){
            chrNum = 23L;
        }
        else if(dummy[1]=="Y"){
            chrNum = 24L;
        }
        else{
            chrNum = Long.parseLong(dummy[1]);
        }
                
        return chrNum;
    }
    
    public StringBuffer getSequence() {
        try{
            BufferedReader chr_reader = new BufferedReader(new InputStreamReader(new FileInputStream(this.getFilePath()+".fa")));
            System.out.println(chr_reader.readLine());
            StringBuffer sb = new StringBuffer(chr_reader.readLine());
            System.out.println(sb.subSequence(10000, 10100));
            this.setSequence(sb);
        }catch(IOException e){
            
        }
        return seq;
    }

    public String getFilePath() {
        return ref.getPath()+File.separator+name;
    }

    public void writeToFile(String fa) throws FileNotFoundException {
        String path = getFilePath();
        writeToPath(path, fa);
    }

    public void setSequence(StringBuffer sb) {
        this.seq = sb;
    }
    
}
