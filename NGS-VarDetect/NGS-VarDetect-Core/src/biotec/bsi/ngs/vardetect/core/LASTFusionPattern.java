/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

/**
 *
 * @author worawich
 */
public class LASTFusionPattern {
    int frontID;
    int backID;
    long overlapBase = 0;
    
    public LASTFusionPattern(){
        
    }

    public void setFrontID(int frontID) {
        this.frontID = frontID;
    }

    public void setBackID(int backID) {
        this.backID = backID;
    }

    public void setOverlapBase(long overlapBase) {
        this.overlapBase = overlapBase;
    }

    public int getFrontID() {
        return frontID;
    }

    public int getBackID() {
        return backID;
    }

    public long getOverlapBase() {
        return overlapBase;
    }
    
} 