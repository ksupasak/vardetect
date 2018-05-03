/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

import java.util.ArrayList;
import java.util.Comparator;

/**
 *
 * @author worawich
 */
public class SVGroupPair{
    
    SVGroup frontSVGroup;
    SVGroup backSVGroup;
    int insertionSize;
    int insertionJunction;
    
    public SVGroupPair(){

    }

    public SVGroup getFrontSVGroup() {
        return frontSVGroup;
    }

    public void setFrontSVGroup(SVGroup frontSVGroup) {
        this.frontSVGroup = frontSVGroup;
    }

    public SVGroup getBackSVGroup() {
        return backSVGroup;
    }

    public void setBackSVGroup(SVGroup backSVGroup) {
        this.backSVGroup = backSVGroup;
        int bpB_F = frontSVGroup.getRPB();
        int bpF_F = frontSVGroup.getRPF();
        int bpF_B = backSVGroup.getRPF();
        int bpB_B = backSVGroup.getRPB();
        this.insertionSize = Math.abs(bpF_B - bpB_F);               // The size of Insert DNA portion
        this.insertionJunction = Math.abs(bpB_B - bpF_F);           // A junction size, ideally it should be one if it is clean insertion
    }

    public int getInsertionSize() {
        // The size of Insert DNA portion
        return insertionSize;
    }

    public int getInsertionJunction() {
        // A junction size, ideally it should be one if it is clean insertion
        return insertionJunction;
    }

    public void setInsertionSize(int insertionSize) {
        this.insertionSize = insertionSize;
    }
    
    public static Comparator<SVGroupPair> InsertSizeComparatorLowToHigh = new Comparator<SVGroupPair>() {

	public int compare(SVGroupPair s1, SVGroupPair s2) {
            int s1Coverage = s1.getInsertionSize();
            int s2Coverage = s2.getInsertionSize();
	   //ascending order (low to high)
            return s1Coverage-s2Coverage;

	   //descending order   (high to low)
//            return s2Coverage-s1Coverage;
        }
    };
    
    public static Comparator<SVGroupPair> InsertJunctionComparatorLowToHigh = new Comparator<SVGroupPair>() {

	public int compare(SVGroupPair s1, SVGroupPair s2) {
            int s1Coverage = s1.getInsertionJunction();
            int s2Coverage = s2.getInsertionJunction();
	   //ascending order (low to high)
            return s1Coverage-s2Coverage;

	   //descending order   (high to low)
//            return s2Coverage-s1Coverage;
        }
    };
    
    public String ShortSummary(){
        
        return null;
    }
}
