package analyzer.variantInfo;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import analyzer.RefSeq.PullRegionsFromRef;
import analyzer.RefSeq.RefSeqParser;
import analyzer.fileWriters.VCFWriter;
import analyzer.transcriptInfo.CDS;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by mwads on 1/26/16.
 */
public class Variant {

    private String Chr;
    private Integer Pos;
    private String Ref;
    private String Alt;
    private String homhet;
    private String spliceInfo;
    private ArrayList<String> Annotations;
    private ArrayList<CDS> transcripts;
    private String GeneName;
    private ArrayList<Double> OriginalMesScores;
    private ArrayList<Double> VariantMesScores;
    private ArrayList<Double> percentDiffList;



    public Variant(String chr, Integer pos, String spliceInfo, String ref, String alt, String homhet){
        this.Chr = chr;
        this.Pos = pos;
        this.Ref = ref;
        this.Alt = alt;
        this.homhet = homhet;
        this.spliceInfo = spliceInfo;
        this.Annotations = new ArrayList<>();
        this.transcripts = new ArrayList<>();


    }

   // Getters and Setters

    public String getChr() {
        return Chr;
    }

    public void setChr(String chr) {
        Chr = chr;
    }

    public Integer getPos() {
        return Pos;
    }

    public void setPos(Integer pos) {
        Pos = pos;
    }

    public String getRef() {
        return Ref;
    }

    public void setRef(String ref) {
        Ref = ref;
    }

    public String getAlt() {
        return Alt;
    }

    public void setAlt(String alt) {
        Alt = alt;
    }

    public String getHomhet() {
        return homhet;
    }

    public void setHomhet(String homhet) {
        this.homhet = homhet;
    }

    public String getSpliceInfo() {
        return spliceInfo;
    }

    public void setSpliceInfo(String spliceInfo) {
        this.spliceInfo = spliceInfo;
    }

    public ArrayList<Double> getOriginalMesScores() {
        return OriginalMesScores;
    }

    public void setOriginalMesScores(ArrayList<Double> originalMesScores) {
        OriginalMesScores = originalMesScores;
    }

    public ArrayList<Double> getVariantMesScores() {
        return VariantMesScores;
    }

    public void setVariantMesScores(ArrayList<Double> variantMesScores) {
        VariantMesScores = variantMesScores;
    }

    // Creating variants to write

    public VariantContextBuilder creatVariantContext(){
        ArrayList<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(this.Ref,true));
        alleles.add(Allele.create(this.Alt,false));
        VariantContextBuilder builder = new VariantContextBuilder("MES Filtering", this.Chr, Long.valueOf(this.Pos), Long.valueOf(this.Pos), alleles);

        builder.attribute("Gene",this.GeneName);

        builder.noGenotypes();
        return builder;
    }

    // Gathering info for variant

    public void addAnnotation(String annotation){
        this.Annotations.add(annotation);
    }

    @Override
    public String toString() {
        return "Variant{" +
                "Chr='" + Chr + '\'' +
                ", Pos=" + Pos +
                ", Ref='" + Ref + '\'' +
                ", Alt='" + Alt + '\'' +
                ", homhet='" + homhet + '\'' +
                ", spliceInfo='" + spliceInfo + '\'' +
                ", Annotations=" + String.valueOf(Annotations.size()) +
                ", transcripts=" + String.valueOf(transcripts.size()) +
                ", GeneName='" + GeneName + '\'' +
                '}';
    }

    public void parseSpliceInfo(RefSeqParser rsp, PullRegionsFromRef prfr){
        String[] geneList = this.spliceInfo.split("\\),");
        for(String gene : geneList) {
            String[] gene2trans = gene.split("\\(");
            this.GeneName = gene2trans[0];
            if (gene2trans.length == 1)
                return;
            System.out.println(this.GeneName + " " + gene2trans[1]);

            String[] transInfo = gene2trans[1].replace(")","").split(",");
            for (String trans : transInfo) {
                String[] info = trans.split(":");
                CDS cds = new CDS(info[0], gene2trans[0]);
                cds.setcDot(info[2]);
                cds.setCDotList(parseCDot(info[2]));
                cds.setExon(info[1]);
                cds.extractCDS(rsp, prfr);
                this.transcripts.add(cds);
            }
        }
    }

    private ArrayList<String> parseCDot(String cDot){
        Pattern p = Pattern.compile("c\\.(\\d+)([-+])(\\d+)([ATGC])>([ACTG])");
        Matcher m = p.matcher(cDot);
        ArrayList<String> info = new ArrayList<>();
        if(m.find()){
            info.add(m.group(1));
            info.add(m.group(2));
            info.add(m.group(3));
            info.add(m.group(4));
            info.add(m.group(5));
        }
        return info;
    }

    public String makeMESSequenceFile(String outFolder, VCFWriter vw) throws Exception{


        VariantContextBuilder vcb = creatVariantContext();

        String threePrime = null;
        String fivePrime = null;

        FileWriter threePrimeFile = null;
        FileWriter fivePrimeFile = null;

        try {
            threePrimeFile = new FileWriter(outFolder + "threePrime.txt");
            fivePrimeFile = new FileWriter(outFolder + "fivePrime.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }

        LinkedList<String> filteredNames = new LinkedList<>();

        for (CDS cds : transcripts){
//            System.out.println("In cds for loop");


            List<String> CDotList = cds.getCDotList();
            if(CDotList.get(1).equals("-")){
//                System.out.println("In -");
                if ( Integer.valueOf(CDotList.get(2)) <= 20) {
//                    System.out.println("In region <= 20: "+CDotList.get(2));
                    threePrime = new String(outFolder+"threePrime.txt");
                    write3Prime(threePrimeFile,cds);
                }
                else{
                    filteredNames.add(cds.getTransName());
                }
            }
            else{
//                System.out.println("In +");
                if(Integer.valueOf(CDotList.get(2)) <= 6) {
//                    System.out.println("In region <= 6: "+CDotList.get(2));
                    fivePrime = new String(outFolder+"fivePrime.txt");
                    write5Prime(fivePrimeFile,cds);
                }
                else{
                    filteredNames.add(cds.getTransName());
                }
            }
        }


        if(filteredNames.size()!=0) {
            StringBuilder sb = new StringBuilder();
            for (String t: filteredNames){
                sb.append(t+",");
            }
            sb.replace(sb.length()-1,sb.length(),"");
            vcb.attribute("Transcripts", sb.toString());
            vw.writeVar(vcb.make());
        }

        if(threePrime!= null && fivePrime != null)
            throw new Exception("The variant shouldn't end up on both ends of the intron.");

        if(threePrime == null) {
            if(fivePrime != null) {
                try {
                    fivePrimeFile.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return fivePrime;
            }
            else{
                return null;
            }
        }
        else{
            try {
                threePrimeFile.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return threePrime;

        }

    }

    private void write3Prime(FileWriter threePrime, CDS cds){
        try {
            List<String> CDotList = cds.getCDotList();
            String originalSeq = cds.getMES3Prime(Integer.valueOf(CDotList.get(0)));
            StringBuilder sb = new StringBuilder(originalSeq);
//            System.out.println("cdot ref="+CDotList.get(3)+" intron at "+CDotList.get(2)+"="+sb.charAt(20-Integer.valueOf(CDotList.get(2))));
            sb.setCharAt(20-Integer.valueOf(CDotList.get(2)), CDotList.get(4).charAt(0));
            System.out.println("original: "+originalSeq+"\nnewseq: "+sb.toString());
            threePrime.write(">original\n" + originalSeq + "\n>newSeq\n" + sb.toString() + "\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void write5Prime(FileWriter fivePrime, CDS cds){
        try {
            List<String> CDotList = cds.getCDotList();
            String originalSeq = cds.getMES5Prime(Integer.valueOf(CDotList.get(0)));
            StringBuilder sb = new StringBuilder(originalSeq);
//            System.out.println("length="+String.valueOf(originalSeq.length())+" position="+String.valueOf(2+CDotList.get(2)));
//            System.out.println("cdot ref="+CDotList.get(3)+" intron at "+CDotList.get(2)+"="+sb.charAt(Integer.valueOf(CDotList.get(2))));
            sb.setCharAt(Integer.valueOf(2 + Integer.valueOf(CDotList.get(2))), CDotList.get(4).charAt(0));
            System.out.println("original: "+originalSeq+"\nnewseq: "+sb.toString());
            fivePrime.write(">original\n" + originalSeq + "\n>newSeq\n" + sb.toString() + "\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void makeModifiedProtein(){
        for(CDS cds : transcripts){
            cds.makeModifiedProtein();
        }

    }

    public int checkMesSignificance(){
        this.percentDiffList = new ArrayList<>();
        int sigCount = 0;
        int likelySigCount = 0;
        int notSigCount = 0;
        for(int i=0; i < OriginalMesScores.size(); i++){
            Double percentDiff = ((Double.valueOf(VariantMesScores.get(i))-Double.valueOf(OriginalMesScores.get(i))) / Double.valueOf(OriginalMesScores.get(i)) * 100);
            percentDiffList.add(percentDiff);
            System.out.println(percentDiff);
            if(percentDiff < -75)
                sigCount++;
            else if(percentDiff < -50)
                likelySigCount++;
            else
                notSigCount++;
        }


        for(String s : this.Annotations){
            System.out.println(s);
        }
        if(sigCount > (likelySigCount+notSigCount))
            if(this.Annotations.get(0).equals("NA")&&this.Annotations.get(1).equals("NA")&&this.Annotations.get(2).equals("NA")) //This checks to make sure that the variant doesn't have any annotations from 1KG, Exac, or Gerp++
                return 2;
            else
                return 1;
        if(sigCount == (likelySigCount+notSigCount) || likelySigCount > (sigCount+notSigCount))
            return 1;
        return 0;
    }
}
