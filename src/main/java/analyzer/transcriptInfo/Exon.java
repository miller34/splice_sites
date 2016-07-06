package analyzer.transcriptInfo;

import analyzer.Utilities.Utilities;

/**
 * Created by mwads on 1/26/16.
 */
public class Exon {

    private StringBuilder Seq;
    private String Start;
    private String end;
    private String Strand;

    public Exon(String Start, String End, StringBuilder Seq, String Strand){
        this.Seq = new StringBuilder(Seq);
        this.Start = Start;
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
        return Start;
    }

    public void setStart(String start) {
        Start = start;
    }

    public String getEnd() {
        return end;
    }

    public void setEnd(String end) {
        this.end = end;
    }


    public Integer getLength(){
        return Seq.length();
    }

}
