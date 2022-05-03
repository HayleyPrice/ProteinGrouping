package utils;

import group.Group;
import group.GroupArray;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import pep.PepArray;
import pep.Peptide;
import prot.ProtArray;
import prot.Protein;
import org.apache.commons.math3.stat.inference.TTest;

/**
 * Provides methods to save output.
 */
public class Save {

    public Save() {

    }
    /**
     *
     * @param fname
     * @param num
     * @param peptides
     */
    public void savePeps(String fname, int num, PepArray peptides) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(fname), false);
            outFile.println("PepName ,ProtNo ,pepType ,Quants ,,,,,,,,,,,,Proteins");
            for (int i = 0; i < peptides.getSize(); i++) {
                Peptide p = peptides.getPep(i);
                outFile.print(p.getPepName() + "," + p.getProtNo()
                        + "," + p.pepType() + ",");
                for (int j = 0; j < num; j++) {
                    outFile.print(p.getQuantVals(j) + ",");
                }
                List<Protein> prots = p.getProtList();
                for (int k = 0; k < prots.size(); k++) {
                    Protein pr = prots.get(k);
                    outFile.print(pr.getProtName() + ",");
                }
                outFile.println();
            }
            outFile.close();
        } catch (Exception e) {System.out.println("Unable to save to " + fname);}
    }
    public void saveMSstatsConf(String op, String fname, int num, int conds, PepArray peptides, GroupArray groups) {
        char condition = '\0';
        int condRuns = num / conds;
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(op + fname), false);
            outFile.println("ProteinName ,PeptideSequence ,PrecursorCharge ,FragmentIon ,"
                    + "ProductCharge ,IsotopeLabelType ,Condition ,BioReplicate , "
                    + "Run ,Intensity");
            for (int i = 0; i < groups.getSize(); i++) {
                Group g = groups.getGroup(i);
                Protein gh = g.getGroupHead();
                List<Protein> prots = g.getProtGroupList();
                String gpNames = "";
                if (!gh.isDiscarded) {                    
                    for (Protein pr : prots) {
                        gpNames += pr.getProtName() + ";";                        
                    }
                    List<Peptide> pepList = gh.getPepList();
                    
                    int qPeps = 0;
                    for (Peptide pep : pepList) {
                        //if (pep.forQuant) {
                        if (pep.isUnique) {
                            qPeps++;
                        }
                    }
                    //if (qPeps >= 1 ) {
                    if (qPeps >= 2 ) {
                        for (Peptide pep : pepList) {
                            for (int run = 0; run < num; run++) {
                                if (run < condRuns) {
                                    condition = 'A';
                                }
                                if (run < 2 * condRuns && run >= condRuns) {
                                    condition = 'B';
                                }
                                if (run < 3 * condRuns && run >= 2 * condRuns) {
                                    condition = 'C';
                                }
                                if (run < 4 * condRuns && run >= 3 * condRuns) {
                                    condition = 'D';
                                }
                                if (pep.isConflicted) {
                                //outFile.print(gpNames + "," + pep.getPepIdent() + "_" + gpNames           //to print all proteins in group
                                outFile.print(gh.getProtName() + "," + pep.getPepIdent() + "_" + gpNames    //to print only head protein for benchmarking
                                        + "," + pep.getCharge() + "," + "NA,NA,L" 
                                        + "," + condition + ",1," + run + "," + 
                                        + pep.getAbundShare(gh.getProtName(), run, num));
                                outFile.println();
                                }
                                else {                        
                                    //outFile.print(gpNames + "," + pep.getPepIdent() + "_" + gpNames           //to print all proteins in group
                                outFile.print(gh.getProtName() + "," + pep.getPepIdent() + "_" + gpNames    //to print only head protein for benchmarking 
                                            + "," + pep.getCharge() + "," + "NA,NA,L" 
                                            + "," + condition + ",1," + run + "," + 
                                            + pep.getQuantVals(run));
                                    outFile.println();                            
                                }
                            }
                        }
                    }
                }
            }
            outFile.close();
        } catch (Exception e) {System.out.println("Unable to save to " + fname);}
    }
    public void saveMSstatsNC(String op, String fname, int num, int conds, PepArray peptides, GroupArray groups) {
        char condition = '\0';
        int condRuns = num / conds;
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(op + fname), false);
            outFile.println("ProteinName ,PeptideSequence ,PrecursorCharge ,FragmentIon ,"
                    + "ProductCharge ,IsotopeLabelType ,Condition ,BioReplicate , "
                    + "Run ,Intensity");
            for (int i = 0; i < groups.getSize(); i++) {
                Group g = groups.getGroup(i);
                Protein gh = g.getGroupHead();
                List<Protein> prots = g.getProtGroupList();
                String gpNames = "";
                if (!gh.isDiscarded) {                    
                    for (Protein pr : prots) {
                        gpNames += pr.getProtName() + ";";                        
                    }
                    List<Peptide> pepList = gh.getPepList();
                    
                    int qPeps = 0;
                    for (Peptide pep : pepList) {
                        //if (pep.forQuant) {
                        if (pep.isUnique) {
                            qPeps++;
                        }
                    }
                    //if (qPeps >= 1 ) {
                    if (qPeps >= 2 ) {
                        for (Peptide pep : pepList) {
                            for (int run = 0; run < num; run++) {
                                if (run < condRuns) {
                                    condition = 'A';
                                }
                                if (run < 2 * condRuns && run >= condRuns) {
                                    condition = 'B';
                                }
                                if (run < 3 * condRuns && run >= 2 * condRuns) {
                                    condition = 'C';
                                }
                                if (run < 4 * condRuns && run >= 3 * condRuns) {
                                    condition = 'D';
                                }
                                if (!pep.isConflicted) {
                                    outFile.print(gh.getProtName() + "," + pep.getPepIdent() + "_" + gpNames    //to print only head protein for benchmarking 
                                            + "," + pep.getCharge() + "," + "NA,NA,L" 
                                            + "," + condition + ",1," + run + "," + 
                                            + pep.getQuantVals(run));
                                    outFile.println();                            
                                }
                            }
                        }
                    }
                }
            }
            outFile.close();
        } catch (Exception e) {System.out.println("Unable to save to " + fname);}
    }
    /**
     *
     * @param fname
     * @param proteins
     */
    public void saveProts(String fname, ProtArray proteins) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(fname), false);
            outFile.println("Protein ,PepNo ,ProtType ,Discarded? ,Peptides");
            for (int i = 0; i < proteins.getSize(); i++) {
                Protein p = proteins.getProt(i);
                outFile.print(p.getProtName() + "," + p.getPepNo() + "," 
                        + p.protType() + "," + p.isDiscarded + ",");      
                List<Peptide> peps = p.getPepList();
                for (Peptide ps : peps) {
                    outFile.print(ps.getPepName() + "," + ps.getProtNo() + ",");
                }
                outFile.println();
            }
            outFile.close();
        } catch (Exception e) {System.out.println("Unable to save to " + fname);}
    }
    public void saveQProt(String op, String fname, GroupArray groups, 
            String quantMethod) {
        for (int i = 0; i < 3; i++) {   //PXD001385/PXD001819
        //for (int i = 0; i < ; i++) { //PXD002099
            try {
                PrintWriter outFile = new PrintWriter(new FileWriter(
                        op + fname + i), false);
                outFile.println("Protein \t0\t0\t0\t1\t1\t1");
                //outFile.println("Protein ");
                //outFile.println("Protein ");
                for (int j = 0; j < groups.getSize(); j++) {
                    Group g = groups.getGroup(j);
                    Protein gh = g.getGroupHead();
                    if (!gh.isDiscarded) {
                        String groupHeadName = gh.getProtName();
                        //outFile.print(groupHeadName);
                        List<Protein> prots = g.getProtGroupList();
                        List<Peptide> pepList = gh.getPepList();
                        int qPeps = 0;
                        for (Peptide pep : pepList) {
                            //if (pep.forQuant) {
                            if (pep.isUnique) {
                                qPeps++;
                            }
                        }
                        //if (qPeps >= 1 && !groupHeadName.isEmpty()) {
                        if (qPeps >= 2 && !groupHeadName.isEmpty()) {
                            outFile.print(groupHeadName);
//                            for (Protein pr : prots) {
//                                String protName = pr.getProtName();
//                                if (!protName.equals(groupHeadName)) {
//                                    outFile.print(";" + protName);
//                                }                        
//                            }
                            outFile.print("\t");
                            if (i == 0) {
                                for (int run = 0; run < 3; run++) {
                                    outFile.print(gh.getQuant(run, quantMethod) + "\t");
                                }
                            }
                            if (i == 1) {
                                for (int run = 3; run < 6; run++) {
                                    outFile.print(gh.getQuant(run, quantMethod) + "\t");
                                }
                            }
                            if (i == 2) {
                                for (int run = 6; run < 9; run++) {
                                    outFile.print(gh.getQuant(run, quantMethod) + "\t");
                                }
                            }
//                            if (i == 3) {
//                                for (int run = 9; run < 12; run++) {
//                                    outFile.print(gh.getQuant(run, quantMethod) + "\t");
//                                }
//                            }
//                            if (i == 4) {
//                                for (int run = 12; run < 15; run++) {
//                                    outFile.print(gh.getQuant(run, quantMethod) + "\t");
//                                }
//                            }
//                            for (int run = 15; run < 17; run++) {
//                                    outFile.print(gh.getQuant(run, quantMethod) + "\t");
//                            }
//                            for (int run = 17; run < 18; run++) {
//                                    outFile.print(gh.getQuant(run, quantMethod));
//                            }
                            for (int run = 9; run < 11; run++) {
                                    outFile.print(gh.getQuant(run, quantMethod) + "\t");
                            }
                            for (int run = 11; run < 12; run++) {
                                    outFile.print(gh.getQuant(run, quantMethod));
                            }
                            outFile.println();
                        }
                    }
                }
                outFile.close();
            } catch (Exception e) {System.out.println("Unable to save to " + 
                    "QPROT\\" + fname + "_" + i);}
        }
    }
    /**
     *
     * @param fname     * @param num
     * @param groups
     */    
    public void saveGroups(String fname, int num, GroupArray groups, 
            String quantMethod) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(fname), false);
            outFile.println("Proteins ,Peptides ,PepsForQuant ,Spike ,BG");
            for (int i = 0; i < groups.getSize(); i++) {
                Group g = groups.getGroup(i);
                Protein gh = g.getGroupHead();
                if ((!gh.isDiscarded) && (gh.isSpike || gh.isBackground)) {
                    String groupHeadName = gh.getProtName();
                    
                    List<Peptide> pepList = gh.getPepList();
                    
                    int qPeps = 0;
                    for (Peptide pep : pepList) {
                        //if (pep.forQuant) {
                        if (pep.isUnique) {
                            qPeps++;
                        }
                    }
                    //if (qPeps >= 1 ) {
                    if (qPeps >= 2 ) {
                        outFile.print(groupHeadName);
                        List<Protein> prots = g.getProtGroupList();
                        for (Protein pr : prots) {
                            String protName = pr.getProtName();
                            if (!protName.equals(groupHeadName)) {
                                outFile.print(";" + protName);
                            }                        
                        }
                        outFile.print(",");
                        outFile.print(pepList.size() + ",");                        
                        outFile.print(qPeps + ",");
                        if (gh.isSpike) {
                                outFile.print(1 + ",0,");
                        }
                        else if (gh.isBackground) {
                                outFile.print("0," + 1 + ",");
                        }
                        else {
                                outFile.print("0,0,");
                        }
                        for (int run = 0; run < num; run++) {
                            outFile.print(gh.getQuant(run, quantMethod) + ",");
                        }
                        outFile.println();
                    }
                }
            }
            outFile.close();
        } catch (Exception e) {
            System.out.println("Unable to save to " + fname);
        }
    }
    public void saveTtest(String op, String fname, int num, GroupArray groups, 
            String quantMethod, int comp) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(
                    op + fname + "_" + comp + ".csv"), false);
            outFile.println("Protein , pVal ,spikeTot ,bgTot ,FDR ,"
                    + "Sensitivity ,qVal");
            //System.out.println(op + "Ttest\\" + fname + "_" + comp + ".csv");
            for (int j = 0; j < groups.getSize(); j++) {
                Group g = groups.getGroup(j);
                Protein gh = g.getGroupHead();
                double p1 = gh.getPval(comp);
                double p2 = gh.getPval(comp);
                if (p1 != p2) {
                }
                else {
                    if (!gh.isDiscarded) {
                        String groupHeadName = gh.getProtName();
                        outFile.print(groupHeadName);
                        List<Protein> prots = g.getProtGroupList();
                        for (Protein pr : prots) {
                            String protName = pr.getProtName();
                            if (!protName.equals(groupHeadName)) {
                                outFile.print(";" + protName);
                            }                        
                        }
                        outFile.print("," + gh.getPval(comp) + ",");
                        if (gh.isSpike) {
                            outFile.print(1 + ",0,");
                        }
                        else if (gh.isBackground) {
                            outFile.print("0," + 1 + ",");
                        }
                        else {
                            outFile.print(",,");
                        }
                        outFile.print( gh.getFDRs(comp) 
                                + "," + gh.getSensitivities(comp) + "," + gh.getqVals(comp));
                        outFile.println();
                    }
                }
            }
            outFile.close();
        } catch (Exception e) {
            System.out.println("Unable to save to " + op + "Ttest\\" 
                    + fname + "_" + comp  + ".csv");
        }
    }
    public void saveQlone(String op, String fname, int num, int conds, 
            PepArray peptides, GroupArray groups, String quantMethod) {
        
        char condition = '\0';
        int condRuns = num / conds;
        for (int c = 0; c < 3; c++) {
            try {
                PrintWriter outFile = new PrintWriter(new FileWriter(op + fname), false);

                for (int i = 0; i < groups.getSize(); i++) {
                    Group g = groups.getGroup(i);
                    Protein gh = g.getGroupHead();
                    String ghName = gh.getProtName();
                    List<Protein> prots = g.getProtGroupList();
                    String gpNames = "";
                    if (!gh.isDiscarded) {                    
                        for (Protein pr : prots) {
                            gpNames += pr.getProtName() + ";";                        
                        }
                        List<Peptide> pepList = gh.getPepList();
                        int qPeps = 0;
                        for (Peptide pep : pepList) {
                            //if (pep.forQuant) {
                            if (pep.isUnique) {
                                qPeps++;
                            }
                        }
                        //if (qPeps >= 1 ) {
                        if (qPeps >= 2 ) {
                            //outFile.print(gpNames + "\t");
                            if (c == 0) {
                                for (int run = 0; run < 3; run++) {
                                    outFile.println(ghName + "\t0\t " + gh.getQuant(run, quantMethod));
                                }
                            }
                            if (c == 1) {
                                for (int run = 3; run < 6; run++) {
                                    outFile.println(ghName + "\t0\t " + gh.getQuant(run, quantMethod));
                                }
                            }
                            if (c == 2) {
                                for (int run = 6; run < 9; run++) {
                                    outFile.println(ghName + "\t0\t " + gh.getQuant(run, quantMethod));
                                }
                            }
                            for (int run = 9; run < 12; run++) {
                                outFile.println(ghName + "\t1\t " + gh.getQuant(run, quantMethod));
                            }
                        }
                    }
                }
            outFile.close();
            } catch (Exception e) {System.out.println("Unable to save to " + fname);}
        }
    }
}

