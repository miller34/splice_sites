package analyzer.transcriptInfo;

import analyzer.RefSeq.PullRegionsFromRef;
import analyzer.RefSeq.RefSeqParser;
import analyzer.Utilities.Utilities;

import java.util.*;

/**
 * This class holds the information for the CDS region in question, extracts the CDS sequences, and pulls the sequences
 * used in MaxEntScan.
 *
 * Created by mwads on 1/26/16.
 */
public class CDS {

    private ArrayList<Exon> Exons;
    private ArrayList<Intron> Introns;
    private String Protein;
    private String geneName;
    private String transName;
    private String cdsStart; // column 6
    private String cdsEnd; // column 7
    private String chr; // column 2 refseq
    private String strand; // column 3
    private String[] exonStarts; // column 9
    private String[] exonEnds; // column 10
    private String[] exonFrames; // column 15
    private StringBuilder seq;
    private String cDot;
    private ArrayList<String> CDotList;
    private String exon;
    private Integer prime; //this will be either 3 or 5 depending on which side of the intron it is
    private Integer exonSpliceMissed; // the number of the exon that is closest to the offending variant
    private String modifiedProtein;



    public CDS(String transName, String gene) {
        this.geneName = gene;
        this.transName = transName;
        this.Exons = new ArrayList<>();
        this.Introns = new ArrayList<>();
        this.seq = new StringBuilder();


    }

    //Getters and Setters
    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand)
    {
        this.strand = strand;
    }

    public String getcDot() {
        return cDot;
    }

    public void setcDot(String cDot) {
        this.cDot = cDot;
    }

    public ArrayList<String> getCDotList() {
        return CDotList;
    }

    public void setCDotList(ArrayList<String> CDotList) {
        this.CDotList = CDotList;
    }

    public String getExon() {
        return exon;
    }

    public void setExon(String exon) {
        this.exon = exon;
    }

    public String getTransName(){
        return this.transName;
    }


    public void extractCDS(RefSeqParser rsp, PullRegionsFromRef prfr) {
        String info = rsp.getRefSeqData(transName);
        String[] data = info.split("\\t");
        this.chr = data[2];
        this.strand = data[3];
        Integer exonNum = Integer.valueOf(data[8]);
        this.cdsStart = data[6];
        this.cdsEnd = data[7];
        this.exonStarts = data[9].split(",");
        this.exonEnds = data[10].split(",");
        this.exonFrames = data[15].split(",");
        this.seq = new StringBuilder();

        boolean firstExon = true;

        if(this.strand.equals("+")){
            extractCDSRegionPosStrand(prfr,exonNum,firstExon);
        }
        else{
            extractCDSRegionNegStrand(prfr,exonNum,firstExon);
        }

    }

    private void extractCDSRegionNegStrand(PullRegionsFromRef prfr, Integer exonNum, boolean firstExon){
        for(int i=exonNum-1;i >= 0;i--){
            if(Integer.valueOf(this.exonFrames[i]) != -1) {
                if(!firstExon && Integer.valueOf(this.cdsStart) < Integer.valueOf(this.exonStarts[i])){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.exonEnds[i],this.exonStarts[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
//                    System.out.println(e.getSeq());
                }
                else if(firstExon){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(cdsEnd), this.strand);
                    Exon e = new Exon(this.cdsEnd,this.exonStarts[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
//                    System.out.println(e.getSeq());

                    firstExon = false;
                }
                else{
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.cdsStart), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.exonEnds[i],this.cdsStart,exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
//                    System.out.println(e.getSeq());
                }

                if(i > 0) {
//                    System.out.println(String.valueOf(Integer.valueOf(this.exonEnds[i - 1]))+":"+String.valueOf(Integer.valueOf(this.exonStarts[i])));
                    StringBuilder intron = prfr.getRegion(this.chr, Integer.valueOf(this.exonEnds[i - 1]), Integer.valueOf(this.exonStarts[i]), this.strand);
//                    System.out.println(intron.toString());
                    Intron I = new Intron(intron, this.exonStarts[i], this.exonEnds[i - 1],this.strand);

                    this.Introns.add(I);
                }
            }
        }

        this.Protein = Utilities.translateProt(this.seq);
//        System.out.println(this.Protein);
    }

    private void extractCDSRegionPosStrand(PullRegionsFromRef prfr, Integer exonNum, boolean firstExon){
        for(int i=0;i < exonNum;i++){
            if(Integer.valueOf(this.exonFrames[i]) != -1) {
                if(!firstExon && Integer.valueOf(this.cdsEnd) > Integer.valueOf(this.exonEnds[i])){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.exonStarts[i],this.exonEnds[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
                }
                else if(firstExon){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(cdsStart), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.cdsStart,this.exonEnds[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());

                    firstExon = false;
                }
                else{
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(this.cdsEnd), this.strand);
                    Exon e = new Exon(this.exonStarts[i],this.cdsEnd,exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
                }

                if(i < exonNum-1) {
                    StringBuilder intron = prfr.getRegion(this.chr, Integer.valueOf(this.exonEnds[i]), Integer.valueOf(this.exonStarts[i + 1]), this.strand);
                    Intron I = new Intron(intron, this.exonEnds[i], this.exonStarts[i + 1],this.strand);

                    this.Introns.add(I);
                }
            }
        }

        this.Protein = Utilities.translateProt(this.seq);
    }

    public String getMES3Prime(Integer pos){
//        if(this.strand.equals("+")){
            StringBuilder sb = new StringBuilder();
            Integer e = -1;
//            if (this.strand.equals("+")){
                e = getPosExon(pos);
//            }
//            else{
//                e = getNegExon(pos);
//            }
//            System.out.println(String.valueOf(e));

            this.exonSpliceMissed = e;
            this.prime = 3;
//            System.out.println(this.Introns.get(e-1).getStart()+":"+this.Introns.get(e-1).getEnd()+"-"+this.Introns.get(e-1).getSeq());

            StringBuilder intron = new StringBuilder(this.Introns.get(e-1).getSeq());
            sb.append(intron.substring(intron.length()-20,intron.length()));
            StringBuilder seq = new StringBuilder(this.Exons.get(e).getSeq());
            sb.append(seq.substring(0,3));
            return sb.toString();
//        }
//        else{
//            StringBuilder sb = new StringBuilder();
//            Integer e = getNegExon(pos);
//            StringBuilder intron = new StringBuilder(this.Introns.get(e-1).getSeq());
//            sb.append(intron.substring(intron.length()-20,intron.length()));
//            StringBuilder seq = new StringBuilder(this.Exons.get(e).getSeq());
//            sb.append(seq.substring(0,3));
//            return sb.toString();
//        }
    }



    public String getMES5Prime(Integer pos){
//        if(this.strand.equals("+")){
            StringBuilder sb = new StringBuilder();
            Integer e = -1;
//            if (this.strand.equals("+")){
                e = getPosExon(pos);
//            }
//            else{
//                e = getNegExon(pos);
//            }
//            System.out.println(e);


            this.exonSpliceMissed = e;
            this.prime = 5;
            StringBuilder seq = new StringBuilder(this.Exons.get(e).getSeq());
            sb.append(seq.substring(seq.length()-3,seq.length()));
            StringBuilder intron = new StringBuilder(this.Introns.get(e).getSeq());
            sb.append(intron.substring(0,6));
            return sb.toString();
//        }
//        else{
//            StringBuilder sb = new StringBuilder();
//            Integer e = getNegExon(pos);
//            StringBuilder seq = new StringBuilder(this.Exons.get(e).getSeq());
//            sb.append(seq.substring(seq.length()-3,seq.length()));
//            StringBuilder intron = new StringBuilder(this.Introns.get(e-1).getSeq());
//            sb.append(intron.substring(0,6));
//            return sb.toString();
//        }
    }

    private Integer getNegExon(Integer Pos){
        Integer total = 0;
        Integer i;
        for(i = 0; i < this.Exons.size(); i++){
//            System.out.println("start: "+this.exonStarts[i]+" end: "+this.exonEnds[i]);
            total += this.Exons.get(i).getLength();
//            System.out.println("Negative: "+String.valueOf(Pos)+"<"+String.valueOf(total));
            if(Pos < total)
                continue;
            break;
        }
        return i;
    }

    private Integer getPosExon(Integer Pos){
        Integer total = 0;
        Integer i;
        for(i = 0; i < this.Exons.size(); i++){
            total += this.Exons.get(i).getLength();
//            System.out.println("Postive: "+String.valueOf(Pos)+">"+String.valueOf(total));
            if(Pos > total)
                continue;
            break;
        }
        return i;
    }


    public String makeModifiedProtein(){
        StringBuilder sb = new StringBuilder();
//        System.out.println(exonSpliceMissed);
//        System.out.println(prime);
        for(int i = 0; i < Exons.size(); i++){
            if(i == exonSpliceMissed && prime == 3){
                sb.append(Introns.get(exonSpliceMissed-1).getSeq());
            }
            sb.append(Exons.get(i).getSeq());
//            System.out.println(Exons.get(i).getSeq());
            if(i == exonSpliceMissed && prime == 5){
                sb.append(Introns.get(exonSpliceMissed).getSeq());
            }
        }
        modifiedProtein = Utilities.translateProt(sb);

        System.out.println(this.Protein);
        System.out.println(this.modifiedProtein);

        return modifiedProtein;
    }


}
