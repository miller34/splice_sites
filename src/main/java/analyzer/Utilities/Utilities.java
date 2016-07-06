package analyzer.Utilities;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

/**
 * Created by mwads on 2/10/16.
 */
public class Utilities {

    public static final String RESET = "\u001B[0m";
    public static final String GREEN = "\u001B[32m";


    public static String getProcessError(Process p){
        InputStream input = p.getErrorStream();
        BufferedReader reader = new BufferedReader(new InputStreamReader(input));
        StringBuilder sb = new StringBuilder();
        String line = null;
        try {
            while( (line = reader.readLine()) != null){
                sb.append(line);
                sb.append(System.getProperty("line.separator"));
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sb.toString();
    }

    public static String getProcessOutput(Process p){
        InputStream input = p.getInputStream();
        BufferedReader reader = new BufferedReader(new InputStreamReader(input));
        StringBuilder sb = new StringBuilder();
        String line = null;
        try {
            while( (line = reader.readLine()) != null){
                sb.append(line);
                System.out.println(line);
                sb.append(System.getProperty("line.separator"));
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sb.toString();
    }

    public static StringBuilder reverseComplement(StringBuilder sequence){
        StringBuilder seq = new StringBuilder(sequence);
        seq = seq.reverse();
        HashMap<String,String> complements = new HashMap<>();
        complements.put("A","T");
        complements.put("T","A");
        complements.put("C","G");
        complements.put("G","C");
        complements.put("N","N");
        for(int j = 0; j < seq.length(); j++) {
            seq.replace(j, j + 1, complements.get(seq.substring(j, j + 1).toUpperCase()));
        }
        return seq;
    }

    public static String translateProt(StringBuilder seq){
        HashMap<String,String> aminoAcids = new HashMap<>();
        aminoAcids = createAminoAcidChart(aminoAcids);
        StringBuilder prot = new StringBuilder();
        for(int i =0; i < seq.length()-2; i+=3){
            String aa = aminoAcids.get(seq.substring(i,i+3).toUpperCase());
            prot.append(aa);
            if(aa.equals("*")){
                break;
            }
        }
        return prot.toString();
    }

    private static HashMap<String,String> createAminoAcidChart(HashMap<String,String> aaChart){
        aaChart.put("TTT","F");
        aaChart.put("TTC","F");
        aaChart.put("TTA","L");
        aaChart.put("TTG","L");
        aaChart.put("TCT","S");
        aaChart.put("TCC","S");
        aaChart.put("TCA","S");
        aaChart.put("TCG","S");
        aaChart.put("TAT","Y");
        aaChart.put("TAC","Y");
        aaChart.put("TAA","*");
        aaChart.put("TAG","*");
        aaChart.put("TGT","C");
        aaChart.put("TGC","C");
        aaChart.put("TGA","*");
        aaChart.put("TGG","W");
        aaChart.put("CTT","L");
        aaChart.put("CTC","L");
        aaChart.put("CTA","L");
        aaChart.put("CTG","L");
        aaChart.put("CCT","P");
        aaChart.put("CCC","P");
        aaChart.put("CCA","P");
        aaChart.put("CCG","P");
        aaChart.put("CAT","H");
        aaChart.put("CAC","H");
        aaChart.put("CAA","Q");
        aaChart.put("CAG","Q");
        aaChart.put("CGT","R");
        aaChart.put("CGC","R");
        aaChart.put("CGA","R");
        aaChart.put("CGG","R");
        aaChart.put("ATT","I");
        aaChart.put("ATC","I");
        aaChart.put("ATA","I");
        aaChart.put("ATG","M");
        aaChart.put("ACT","T");
        aaChart.put("ACC","T");
        aaChart.put("ACA","T");
        aaChart.put("ACG","T");
        aaChart.put("AAT","N");
        aaChart.put("AAC","N");
        aaChart.put("AAA","K");
        aaChart.put("AAG","K");
        aaChart.put("AGT","S");
        aaChart.put("AGC","S");
        aaChart.put("AGA","R");
        aaChart.put("AGG","R");
        aaChart.put("GTT","V");
        aaChart.put("GTC","V");
        aaChart.put("GTA","V");
        aaChart.put("GTG","V");
        aaChart.put("GCT","A");
        aaChart.put("GCC","A");
        aaChart.put("GCA","A");
        aaChart.put("GCG","A");
        aaChart.put("GAT","D");
        aaChart.put("GAC","D");
        aaChart.put("GAA","E");
        aaChart.put("GAG","E");
        aaChart.put("GGT","G");
        aaChart.put("GGC","G");
        aaChart.put("GGA","G");
        aaChart.put("GGG","G");

        return aaChart;
    }
}
