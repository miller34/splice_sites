package analyzer.maxEntScan;

import analyzer.Utilities.Utilities;
import analyzer.fileWriters.VCFWriter;
import analyzer.variantInfo.Variant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;

/**
 * This class runs MaxEntScan to score each splice site which contains a variant.
 *
 * Created by mwads on 3/2/16.
 */
public class MESRunner {

    private String path;
    private Boolean isEmpty = true;

    private String[] results = null;
    private String error;

    public Boolean IsEmpty() {
        return isEmpty;
    }

    /**
     * Creates the sequence file in the variant class and scores the resulting sequences with MaxEntScan and
     * gives the score back to the variant file.
     *
     * @param var
     * @param outFolder
     * @param vw
     * @param algorithmPath
     * @throws Exception
     */
    public MESRunner(Variant var,String outFolder,VCFWriter vw, String algorithmPath) throws Exception{

        this.path = algorithmPath;
        String file = var.makeMESSequenceFile(outFolder, vw);
        if (file == null){
            this.isEmpty = true;
        }
        else if (file.contains("three")){
            run3Prime(file);
            this.isEmpty = false;
            assignScores(var);
        }
        else{
            run5Prime(file);
            this.isEmpty = false;
            assignScores(var);
        }

    }

    /**
     * This function assigns the score list to the variant storage object.
     *
     * @param var
     * @throws Exception
     */
    private void assignScores(Variant var) throws Exception{
        ArrayList<ArrayList<Double>> scores = getScoreLists();

        var.setOriginalMesScores(scores.get(0));
        var.setVariantMesScores(scores.get(1));

    }

    /**
     * Scores the splice site at the 5' end of the intron.
     *
     * @param file
     * @throws Exception
     */
    public void run5Prime(String file) throws Exception{
        System.out.println(Utilities.GREEN+"Running Max Ent Scan 5\'"+ Utilities.RESET);
        try {
            String[] call = new String[]{"perl",this.path+"score5.pl",file};
            ProcessBuilder pb = new ProcessBuilder(call);

            Process p = pb.start();
            results = Utilities.getProcessOutput(p).split("\n");
            error = Utilities.getProcessError(p);

            p.waitFor();
            p.destroyForcibly();

            Files.deleteIfExists(new File(file).toPath());


            if(!error.isEmpty()){
                throw new Exception("MaxEntScan threw the following error: "+error);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }

    /**
     * Scores the 3' end of the intron.
     *
     * @param file
     * @throws Exception
     */
    public void run3Prime(String file) throws Exception{
        System.out.println(Utilities.GREEN+"Running Max Ent Scan 3\'"+ Utilities.RESET);
        try {

            String[] call = new String[]{"perl", this.path+"score3.pl",file};
            ProcessBuilder pb = new ProcessBuilder(call);

            Process p = pb.start();
            results = Utilities.getProcessOutput(p).split("\n");
            error = Utilities.getProcessError(p);

            p.waitFor();
            p.destroyForcibly();

            Files.deleteIfExists(new File(file).toPath());

            if(!error.isEmpty()){
                throw new Exception("MaxEntScan threw the following error: "+error);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }

    /**
     *Gets the string of the list of scores, intended for testing purposes.
     *
     * @return String
     */
    public String getScores(){
        StringBuilder sb = new StringBuilder();
        for(int i=1;i<results.length;i+=2){
            String[] variant = results[i].split("\\s");
            sb.append(variant[1]+",");
        }
        sb.replace(sb.length()-1,sb.length(),"");
        return sb.toString();
    }

    /**
     * Gets the list of MaxEntScan scores.
     *
     * @return
     * @throws Exception
     */
    private ArrayList<ArrayList<Double>> getScoreLists() throws Exception{
        ArrayList<ArrayList<Double>> scores = new ArrayList<>();
        ArrayList<Double> originalList = new ArrayList<>();
        ArrayList<Double> variantList = new ArrayList<>();

        if(results==null)
            throw new Exception("Results is null.");

        for(int i=0;i<results.length;i+=2){
            Double original = Double.valueOf(results[i].split("\\s")[1]);
            Double variant = Double.valueOf(results[i+1].split("\\s")[1]);

            originalList.add(original);
            variantList.add(variant);

        }

        scores.add(originalList);
        scores.add(variantList);

        return scores;
    }
}
