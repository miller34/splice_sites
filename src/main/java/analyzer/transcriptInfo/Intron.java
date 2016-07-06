package analyzer.transcriptInfo;

import analyzer.Utilities.Utilities;

/**
 * Created by mwads on 1/26/16.
 */
public class Intron {

    private StringBuilder Seq;
    private String start;
    private String end;
    private String Strand;


    public Intron(StringBuilder Seq, String Start, String End, String Strand){
//        if(Strand.equals("-")){
//            this.Seq = new StringBuilder(analyzer.Utilities.reverseComplement(Seq));
//        }
//        else {
            this.Seq = new StringBuilder(Seq);
//        }
        this.start = Start;
        this.end = End;
        this.Strand = Strand;
    }

    public String getSeq() {
        return Seq.toString();
    }

    public void setSeq(StringBuilder seq) {
        Seq = seq;
    }

    public String getStart() {
        return start;
    }

    public void setStart(String start) {
        this.start = start;
    }

    public String getEnd() {
        return end;
    }

    public void setEnd(String end) {
        this.end = end;
    }


}
