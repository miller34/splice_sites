package analyzer.RefSeq;

import analyzer.Utilities.Utilities;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Created by mwads on 2/10/16.
 */
public class PullRegionsFromRef {

    private String Path2Ref;
    private String Path2Samtools;

    public PullRegionsFromRef(String Path2Ref, String Path2Samtools){
        StringBuilder sb = new StringBuilder(Path2Ref);
        if (!Path2Ref.endsWith("/")){
            sb.append("/");
        }

        this.Path2Ref = sb.toString();
        this.Path2Samtools = Path2Samtools;
    }

    public StringBuilder getRegion(String chr, Integer start, Integer end, String strand){

        try {

//            System.out.println(String.valueOf(start)+":"+String.valueOf(end));
            Integer startPos = start+1;
            String[] call = null;
                call = new String[]{this.Path2Samtools, "faidx", this.Path2Ref + chr + ".fa", chr + ":" + String.valueOf(startPos) + "-" + String.valueOf(end)};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();

            StringBuilder sb = getFastaSequence(p);
            Utilities.getProcessError(p);

            p.waitFor();

            p.destroyForcibly();
            p.waitFor();

            if(strand.equals("-")){
                sb = Utilities.reverseComplement(sb);
            }

            return sb;

        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;
    }


    private StringBuilder getFastaSequence(Process p){
        BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
        StringBuilder sb = new StringBuilder();
        try {
            String line = null;
            int pos = 0;
            while( (line = reader.readLine()) != null) {
                if (pos != 0) { //The reason for this is to skip the fasta headerline
                    sb.append(line);
                }
                pos += 1;
            }

           reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sb;
    }
}
