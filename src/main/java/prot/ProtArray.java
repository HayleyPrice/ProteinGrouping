package prot;

import pep.Peptide;
import group.Group;
import group.GroupArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Iterator;
import org.apache.commons.math3.stat.inference.TTest;

/**
 * Creates an array of Protein objects.
 */
public class ProtArray {

    private ArrayList<Protein> proteins;

    public ProtArray() {

        proteins = new <Protein>ArrayList();

    }

    /**
     *
     * @param input
     */
    public void buildProtArrayFromProgenesisPepIon(List<String> input) {
        Protein tempProt = null;
        String[] pepProperties = null;
        String splitBy = ",";
        String prot = "";

        for (String line : input) {
            pepProperties = line.split(splitBy);
            prot = pepProperties[10];
            tempProt = checkProt(prot);
        }
    }
    /**
     *
     */
    public void discardNonIdent() {
        Iterator<Protein> iter = proteins.iterator();
        while (iter.hasNext()) {
            Protein p = iter.next();
            if (p.getPepNo() == 0) {
                iter.remove();
            }
        }
    }
    /**
     *
     */
    public void discardSingleIdent() {
        for (Protein p : proteins) {
            if (p.getPepNo() < 1) {       // 1 unique pepetide for identification
            //if (p.getPepNo() < 2) {         // 2 unique pepetides for identification
                p.makeDiscarded();
                List<Peptide> peps = p.getPepList();
            }
        }
    }
    /**
     *
     */
    public void orderProtsByPepCount() {
        // Sort proteins by number of peptides mapped to them
        // Descending
        Collections.sort(proteins,
            (protein1, protein2) -> protein2.getPepNo()
                    - protein1.getPepNo());
        for (Protein p : proteins) {
            List<Peptide> peps = p.getPepList();

            Collections.sort(peps,
                    (pep1, pep2) -> pep1.getProtNo()
                    - pep2.getProtNo());
        }
    }
    /**
     *
     * @param pg
     */
    public void setDistinctProts(GroupArray pg) {
        for (Protein p : proteins) {
            List<Peptide> peps = p.getPepList();
            for (Peptide pep : peps) {
                if (pep.isUnique && !p.isDistinct) {
                    p.makeDistinct();
                    Group gp = new Group(p);
                    pg.addGroup(gp);
                    for (Peptide pp : peps) {
                        pp.makeClaimed();
                        pp.fromDistinct();
                    }
                }
            }
        }
    }
    /**
     *
     * @param pg
     */
    public void setSameSetProts(GroupArray pg) {
        for (Protein p : proteins) {
            List<Peptide> peps = p.getPepList();
            for (Peptide pep : peps) {
                if (!pep.isClaimed && !p.isAssigned) { 
                    p.makeSameSet();
                    Group gp = new Group(p);
                    pg.addGroup(gp);
                    for (Peptide pp : peps) {
                        pp.makeSameSet();
                        pp.makeClaimed();
                        if (!pp.fromDistinct) {
                            pp.makeResolved();
                        }
                    }
                }
                if (pep.fromSameSet && !p.isAssigned) {
                    //Find the group head
                    Group gp = pg.getHeadPeps(pep);
                    if (gp.samePeps(peps)) {
                        gp.addProtToGroup(p);
                    }
                    p.makeSubSet();
                    for (Peptide pp : peps) {
                        pp.makeClaimed();
                        pp.makeSubSet();
                    }
                }
            }
        }
    }
    /**
     *
     * @param pg
     */
    public void setMutSubProts(GroupArray pg) {
        for (Protein p : proteins) {
            List<Peptide> peps = p.getPepList();
            for (Peptide pep : peps) {
                if (pep.isClaimed && !p.isAssigned) {
                    //System.out.println("Pep - " + pep.getPepName());
                    List<Group> gps = new ArrayList<>();
                    for (Peptide pp : peps) {                        
                        pp.makeMutSub();
                        //System.out.println("Pep - " + pp.getPepName());
                        Group g = pg.getHeadPeps(pp);
                        if (g != null) {
                            //System.out.println("ghp - " + pg.getHeadPeps(pp).getGroupHead().getProtName());
                            //System.out.println("groupHead = " + g.getGroupHead().getProtName());
                            gps.add(g);
                            for (Group gp : gps) {
                                gp.addProtToGroup(p);
                            }
                        }
                        else {
                            g = new Group(p);
                        }
                    }
                    p.makeMutSub();
                }
            }
        }
    }
    /**
     *
     * @param pg
     */
    public void discardProts(GroupArray pg) {
        Group gp = null;
        for (Protein p: proteins) {
            if (p.isSameSet) {
                int headPepNo = p.getPepNo();
                gp = pg.getGroup(p);
                if (gp != null) {
                    List<Protein> prots = gp.getProtGroupList();
                    for (Protein pr : prots) {
                        List<Peptide> peps = pr.getPepList();
                        if (headPepNo > pr.getPepNo()) {
                            pr.makeDiscarded();
                            for (Peptide pep : peps) {
                                pep.unClaim();
                                pep.makeResolved();
                            }
                        }
                    }
                }
            }
            if (p.isMutSub) {
                p.makeDiscarded();
                List<Peptide> peps = p.getPepList();
                for (Peptide pep : peps) {
                    if (pep.isClaimed) {
                        pep.makeResolved();
                        pep.unDistinct();
                    }
                }
            }
        }
    }
    /**
     *
     * @param method
     * @param num
     */
    public void setQuants(int num) {        
        for (Protein p: proteins) {
            if (p.isHeadProt && !p.isDiscarded) {
                p.setQuantSum(num);
                p.setQuantHiN(num);
                p.setQuantHi3NonConSameAccrossRuns(num);
            }
        }
    }
    public void setQuantProgenesisHi3NonConflicting(int num) {
        for (Protein p: proteins) {
            p.setQuantProgenesisHi3NonConflicting(num);
        }
    }
    /**
     *
     * @param pt
     * @return
     */
    public Protein retProt(String pt) {
        for (Protein pr : proteins) {
            if (pt.equals(pr.getProtName())) {
                return pr;
            }
        }
        return null;
    }
    /**
     *
     * @return
     */
    public int getSize() {
        return proteins.size();
    }
    /**
     *
     * @param index
     * @return
     */
    public Protein getProt(int index) {
        Protein tempProt = null;
        for (int i = 0; i < proteins.size(); i++) {
            tempProt = proteins.get(index);
        }
        return tempProt;
    }

    private Protein checkProt(String prot) {
        Protein tempProt = null;
        if (!checkProts(prot)) {
            tempProt = newProt(prot);
        }
        else {
            tempProt = retProt(prot);
        }
        return tempProt;
    }
    /**
     *
     * @param prn
     * @return
     */
    private boolean checkProts(String prn){
        for (Protein pr : proteins) {
            if (prn.equals(pr.getProtName())) {
                return true;
            }
        }
        return false;
    }
    private Protein newProt(String pt) {
        Protein tmpProt = new Protein(pt);
        // Add it to the block of proteins
        proteins.add(tmpProt);
        return tmpProt;
    }
}
