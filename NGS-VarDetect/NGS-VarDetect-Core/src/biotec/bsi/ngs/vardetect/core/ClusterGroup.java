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
    
    public ClusterGroup(){
        this.clusterRead = new ArrayList();
    }
    
    public void addShortgunRead(ShortgunSequence readSS){
        this.clusterRead.add(readSS);
    }
    
    public ArrayList<ShortgunSequence> getShortgunRead(){
        return this.clusterRead;
    }
}
