package prot;

import pep.Peptide;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author hayle
 */
public class Protein {

    private String protName;
    private List<Peptide> pepList, NCpepList;
    private List<Double> allPepsQuants, hiNquantsNC, HiNnonConquants, 
            progenesisHiNquantsNC;
    private double[] pVals, FDRs, sensitivities, qVals;
    public boolean isDistinct, isSameSet, isSubset, isMutSub, isAssigned, 
            isHeadProt, isDiscarded, isGroupMember, isSpike, isBackground;
    private int pepNo;
    private String protType;
//    private List<Peptide> resolvedList;
    private Double aveHiNnonConquants;

    
    /**
    * Creates a Protein object.
    *
    * @param    prot  the protein accession.
    */
    public Protein(String prot) {
        
        //PXD001385
        //String background = "HUMAN";    
        //String spike = "ECOLI";
        
        //PXD001819
        //String background = "YEAST";    
        //String spike = "HUMAN";
        
        //PXD002099
        //String background = "YEAST";    
        //String spike = "HUMAN";
        
        String background = "";    
        String spike = "";
        
        this.protName = prot;
        this.pepList = new ArrayList<>();
        this.NCpepList = new ArrayList<>();
        
        this.allPepsQuants = new ArrayList<>();
        this.hiNquantsNC = new ArrayList<>();
        this.HiNnonConquants = new ArrayList<>();
        this.progenesisHiNquantsNC = new ArrayList<>();
        
        this.pVals = new double[2];
        this.FDRs = new double[2];
        this.sensitivities = new double[2];
        this.qVals = new double[2];

        this.pepNo = 0;
        this.aveHiNnonConquants = 0.0;

        isDistinct = false;
        isSameSet = false;
        isMutSub = false;
        isSpike = false;
        isBackground = false;

        isAssigned = false;
        isDiscarded = false;
        isHeadProt = false;
        protType = "";
        
        if (prot.contains(spike)) {
            isSpike = true;
        }
        if (prot.contains(background)) {
            isBackground = true;
        }
    }
    public String getProtName() {
        return this.protName;
    }
    public int getPepNo() {
        return this.pepNo;
    }
    public List getPepList() {
        return this.pepList;
    }
    public void incPepNo() {
        this.pepNo++;
    }
    public void addPepList(Peptide pep) {
        this.pepList.add(pep);        
    }
    public void addNCPepList(Peptide pep) {
        this.NCpepList.add(pep);        
    }    
    public void makeDistinct() {
        this.isDistinct = true;
        this.isAssigned = true;
        this.protType = "Distinct";
    }
    public void makeSameSet() {
        this.isSameSet = true;
        this.isAssigned = true;
        this.protType = "SameSet";
    }
    public void makeSubSet() {
        this.isSubset = true;
        this.isAssigned = true;
        this.protType = "SubSet";
    }
    public void makeMutSub() {
        this.isMutSub = true;
        this.isAssigned = true;
        this.protType = "MutSub";
    }
    public void makeDiscarded() {
        this.isDiscarded = true;
    }
    /**
     * Set quant values using all peptides
     * @param num   Number of runs per condition
     */
    public void setQuantSum(int num) {
        List<Peptide> peps = this.pepList;
        // for each of the runs, i
        for (int run = 0; run < num; run++) {
            Double tempVal = 0.0;
            // peptide quant value for run i for each peptide in list
            // is added to protien quant value i
            for (Peptide pep : peps) {
                // Uncomment for only non-conflicted peptides
                if (!pep.isConflicted) {
                    tempVal = tempVal + pep.getQuantVals(run);
                }
            }
            this.allPepsQuants.add(run, tempVal);
        }        
    }
    /**
     * Sets protein quant values from top 3 most abundant unique/resolved peptides
     * NB maybe different top 3 peptides for each run
     * 
     * @param num   Number of runs per condition
     */
    public void setQuantHiN(int num) {

        List<Peptide> pepL = this.NCpepList;
                
        // Removes any conflicted peptides from list to be used for quantification
        List<Peptide> uniquePeps = removeConflictedPeps(pepL);
        
        for (int i = 0; i < num; i++) {
            List<Double> runPepQuants = new ArrayList<>();
            for (Peptide pep : uniquePeps) {
                runPepQuants.add(pep.getQuantVals(i));
            }
            Collections.sort(runPepQuants, Collections.reverseOrder());

            Double tempVal = 0.0;
            int hiN = 3;
            int vals = runPepQuants.size();
            int loop = hiN;
            if (vals < hiN) {
                loop = vals;
            }
            for (int j = 0; j < loop; j++) {
                tempVal = tempVal + runPepQuants.get(j);                
            }
            this.hiNquantsNC.add(i, tempVal);
        }
    }
    /**
     * Sets protein quant values from top 3 most abundant unique/resolved peptides
     * NB uses same peptides for all runs after ordering the list by the peptide's
     * average abundance across all runs
     * 
     * @param num   Number of runs per condition
     */
    public void setQuantHi3NonConSameAccrossRuns(int num) {
        // Removes any conflicted peptides from list to be used for quantification
        List<Peptide> uniquePeps = this.NCpepList;
        // Orders the peptide list by their average abundance across all runs
        orderPepList(uniquePeps);
        this.HiNnonConquants = getHiNQuants(uniquePeps, num);
        //setAveHiNnonConquants(this.HiNnonConquants);
    }
    private List<Peptide> removeConflictedPeps(List<Peptide> pepList) {
        Iterator<Peptide> iter = pepList.iterator();
        while (iter.hasNext()) {
            Peptide pep = iter.next();
            if (pep.isConflicted) {
                iter.remove();
            }
        }
        return pepList;
    }
    private List<Peptide> orderPepList(List<Peptide> pepList) {
        Collections.sort(pepList,
            (pepetide1, pepetide2) -> pepetide2
                    .getAveAbund().compareTo(pepetide1.getAveAbund()));
        return pepList;
    }
    private List<Double> getHiNQuants(List<Peptide> pepList, int num) {
        List<Double> ncQuants = new ArrayList<>();
        for (int run = 0; run < num; run++) {
            int hiN = 3;
            int pepLNo = pepList.size();
            int loop = hiN;

            if (pepLNo < hiN) {
                loop = pepLNo;
            }
            Double tempVal = 0.0;
            for (int j = 0; j < loop; j++) {         
                tempVal = tempVal + pepList.get(j).getQuantVals(run);
            }
            ncQuants.add(run, tempVal);
        }
        return ncQuants;
    }
    public Double getHiNnonConquants(int run) {
        Double val = 0.0;
        for (Double quant : this.HiNnonConquants) {
            val = this.HiNnonConquants.get(run);
        }
        return val;
    }
    public void setQuantProgenesisHi3NonConflicting(int num) {
        List<Peptide> pepL = orderPepList(this.pepList);
        
        for (int run = 0; run < num; run++) {
            int hiN = 3;
            int pepLNo = pepList.size();
            int loop = hiN;

            if (pepLNo < hiN) {
                loop = pepLNo;
            }
            Double tempVal = 0.0;
            for (int j = 0; j < loop; j++) {
                Peptide pep = pepList.get(j);
                //if (pep.forQuant) {
                    if (!pep.isConflicted) {
                        tempVal = tempVal + pep.getQuantVals(run);
                    }
                    else {                    
                        Double val = pep.getAbundShare(this.protName, run, num);
                        tempVal = tempVal + val;
                    }
                //}
            }
            tempVal = tempVal / hiN;
            this.progenesisHiNquantsNC.add(run, tempVal);
        }
    }
    public String protType() {
        return this.protType;
    }
    public void makeHeadProt() {
        this.isHeadProt = true;
    }
    public double getQuant(int run, String quantMethod) {
        Double quant = 0.0;
        if (quantMethod.equals("HiN")) {
            quant = this.hiNquantsNC.get(run);
        }
        if (quantMethod.equals("sum")) {
            quant = this.allPepsQuants.get(run);
        }
        if (quantMethod.equals("PHN")) {
            quant = this.progenesisHiNquantsNC.get(run);
        }
        return quant;
    }
    public void makeGroupMember() {
        this.isGroupMember = true;
    }
    public void addPvals(int comp, double pval) {        
        this.pVals[comp] = pval;
    }
    public Double getPval(int comp) {
        Double pVal = 0.0;
        for (Double val : this.pVals) {
            pVal = this.pVals[comp];
        }
        return pVal;
    }
    public void addFDRs(int comp, double FDR) {
        this.FDRs[comp] = FDR;
    }
    public Double getFDRs(int comp) {
        Double FDR = 0.0;
        for (Double r : this.FDRs) {
            FDR = this.FDRs[comp];
        }
        return FDR;
    }
        public void addSensitivities(int comp, double s) {
        this.sensitivities[comp] = s;
    }
    public Double getSensitivities(int comp) {
        Double sensitivity = 0.0;
        for (Double s : this.sensitivities) {
            sensitivity = this.sensitivities[comp];
        }
        return sensitivity;
    }
    public void addqVals(int comp, double qVal) {
        this.qVals[comp] = qVal;
    }
    public Double getqVals(int comp) {
        Double qVal = 0.0;
        for (Double q : this.qVals) {
            qVal = this.qVals[comp];
        }
        return qVal;
    }
}
//outFile.println(", Run ,pepprecursor ,Intensity ,ProteinName , "
//                    + "PeptideSequence ,PrecursorCharge ,Condition ,BioReplicate ,"
//                    + "Fragmentation ,ProductCharge ,IsotopeLabelType");