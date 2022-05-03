package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Provides methods to read in data.
 * Methods provided for Progenesis peptide ion and peptide output
 * and mzq file.
 */
public class ReadIn {

    /**
     *
     */
    public ReadIn() {
        
    }
    /**
     * Provides methods to read in Progenesis peptide ion data.
     *
     * @param   fileName    the name and location of the input file
     * @return  data        lines of csv data
     */
    public static List<String> readInCSV(String fileName) {
        int count = 1;
        String tempLine = "";
        List<String> data = new ArrayList<String>();
        try {
            BufferedReader inputFile =
                    new BufferedReader(new FileReader(fileName));
            while ((tempLine = inputFile.readLine()) != null) {
                if (count <= 3 ) {
                    count++;
                }
                else {
                    data.add(tempLine);
                    //System.out.println(tempLine);
                }
            }
        } catch (IOException | NumberFormatException e) {
            System.out.println("Unable to read " + fileName);
        }
    return data;
    } 
}
