/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;

/**
 *
 * @author worawich
 */
public class ClusterGroup {
    
    private ArrayList<ShortgunSequence> clusterRead;
    private ArrayList<String> readNameList;                             // use for local alignment
    
    public ClusterGroup(){
        this.clusterRead = new ArrayList();
        this.readNameList = new ArrayList();
    }
    
    public void addShortgunRead(ShortgunSequence readSS){
        this.clusterRead.add(readSS);
    }
    
    public ArrayList<ShortgunSequence> getShortgunRead(){
        return this.clusterRead;
    }
    
    public void addReadName(String inName){
        this.readNameList.add(inName);
    }
}
