package proteininference;

import group.GroupArray;
import pep.PepArray;
import prot.ProtArray;
import utils.ReadIn;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import utils.Save;

/**
 *
 * @author hayle
 */
public class ProteinInference {

    /**
     *
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        List<String> inputData = new ArrayList<String>();
        //List<String> quantData = new ArrayList<String>();

        PepArray peptides = new PepArray();
        ProtArray proteins = new ProtArray();
        GroupArray protGroups = new GroupArray();
        
         
// **AlterPepArray: 77, Save: 195 for dataset parmeters**
// (All same for Ch2 spike-in data)
        int repNum = 12;    // the number of runs 
        int condNum = 4;    // the number of conditions
        
        //////////////////////////////////////////////////////////////////
        String dataSet = "PXD001385";
        //String fileName = "PXD001385_P_PeptideIon_HiN_NG_23.09.16";
        String fileName = "PXD001385_test";
        
        //String dataSet = "PXD001819";
        //String fileName = "PXD001819_P_PeptideIon_HiN_NG_20.09.16";
        
        //String dataSet = "PXD002099";
        //String fileName = "PXD002099_P_PeptideIon_HiN_NG_25.05.17";
        
        // Set the protein quantification method
        //String quantMethod = "sum";
        
        String quantMethod = "HiN";
        //String quantMethod = "PHN";       
        
        String pepIon = "sumCharge";
        //String pepIon = "sumMods";
        //String pepIon = "sumBoth";
        //String pepIon = "seperate"; 
        
        String progNorm = "_ProgNormPepAb";      // ** Alter line 77 in PepArray **        
        //String progNorm = "_RawPepAb";
        
        //String up = "_1up";                     // ** Alter lines 76, 142, 221, 302 and 411 in save **
        String up = "_2up"; 

        String quantType = "_unique";           // ** Alter lines 71, 137, 216, 297 and 406 in save **
        //String quantType = "_forQuant";
        
        String inputPath = "inputFiles\\";
        String outputPath = "outputFiles\\29.11.20\\";
        
        String outFile = fileName + progNorm + "_" + quantMethod + "_" + pepIon + quantType + up;

        inputData = ReadIn.readInCSV(inputPath + fileName + ".csv");

        peptides.buildPepArrayFromProgenesisPepIon(inputData, repNum, pepIon);
        proteins.buildProtArrayFromProgenesisPepIon(inputData);

        // Maps the proteins and peptides to each other
        peptides.assignProtList(proteins);

        // Discards proteins mapped by 0 peptide
        proteins.discardNonIdent();
        //proteins.discardSingleIdent();

        // Orders the arrays by the number of proteins/peptides mapped
        // Proteins - decreasing peptide number
        // Peptides - increasing protein number
        peptides.orderPepsByProtCount();
        proteins.orderProtsByPepCount();

        // Assigns proteins and peptides
        peptides.setUniquePeptides();
        proteins.setDistinctProts(protGroups);
        proteins.setSameSetProts(protGroups);
        proteins.setMutSubProts(protGroups);
        proteins.discardProts(protGroups);
        peptides.setConflictedPeptides();

        // Quantifies proteins
        proteins.setQuants(repNum);
        
        peptides.setAbundShare(repNum);
        proteins.setQuantProgenesisHi3NonConflicting(repNum);

        //Saves outputs as csv files
        Save saveAs = new Save();
        //saveAs.savePeps(outputPath + "tTest\\" + dataSet + "\\" + "peps_" + outFile + ".csv", repNum, peptides);
        //saveAs.saveProts(outputPath + "tTest\\" + dataSet + "\\" + "prots_" + outFile + ".csv", proteins);
        //saveAs.saveGroups(outputPath + "tTest\\" + dataSet + "\\" + "tTestIn_" + outFile + ".csv", repNum, protGroups, quantMethod);
        saveAs.saveMSstatsConf(outputPath + "MSstats\\", "MSstats_" + outFile + "_R.csv", repNum, condNum, peptides, protGroups);
        //saveAs.saveMSstatsNC(outputPath + "MSstats\\" + dataSet + "\\", "MSstatsIn_" + outFile + ".csv", repNum, condNum, peptides, protGroups);
        //saveAs.saveQProt(outputPath + "QPROT\\" + dataSet + "\\", "QPROTin" + outFile, protGroups, quantMethod);   
   }
}
