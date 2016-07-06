package analyzer;

import analyzer.maxEntScan.MESRunner;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import analyzer.AnnovarRunners.AnnovarRunner;
import analyzer.RefSeq.PullRegionsFromRef;
import analyzer.RefSeq.RefSeqParser;
import analyzer.Utilities.Utilities;
import analyzer.annovarParsers.GeneParser;
import analyzer.annovarParsers.GeneralAnnotationParser;
import analyzer.fileWriters.VCFWriter;
import analyzer.fileWriters.annovarWriter;
import analyzer.variantInfo.Variant;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by mwads on 1/27/16.
 */
public class SpliceRunner {

    private String annovar;
    private String input;
    private String human;
    private String outputFolder;
    private String ref;
    private String refSeq;
    private TreeMap<String,Variant> vars;
    private String SamtoolsPath;
    private String algorithm;
    private String algoritmPath;


    public SpliceRunner(Namespace res, String algorithm){
        this.annovar = res.getString("Annovar");
        this.input = res.getString("Input");
        this.human = res.getString("human");
        this.outputFolder = res.getString("Output");
        if (!this.outputFolder.endsWith("/")){
            StringBuilder sb = new StringBuilder(this.outputFolder+"/");
            this.outputFolder = sb.toString();
        }
        this.ref = res.getString("Ref");
        this.refSeq = res.getString("analyzer/RefSeq");
        this.SamtoolsPath = res.getString("Samtools");
        this.algorithm = algorithm;
        this.algoritmPath = res.getString("AlgorithmPath");

    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("analyzer.SpliceRunner{" +
                "annovar='" + annovar + '\'' +
                ", input='" + input + '\'' +
                ", human='" + human + '\'' +
                ", outputFolder='" + outputFolder + '\'' +
                ", ref='" + ref + '\'' +
                ", refSeq='" + refSeq + "\'vars=");
//        for(Map.Entry<String, Variant> entry : vars.entrySet()){
//            sb.append("\n\tkey= "+entry.getKey()+"\t"+entry.getValue().toString());
//        }

        sb.append(", SamtoolsPath='" + SamtoolsPath + '\'' +
                ", algorithm='" + algorithm + '\'' +
                '}');
        return sb.toString();
    }

    public void run() throws Exception{
        String varfunc = convertAnnovar();
        String newFile = geneBasedAnnotation(varfunc);
        runAnnotations(newFile);
        RefSeqParser rsp = new RefSeqParser(this.refSeq);
        PullRegionsFromRef prfr = new PullRegionsFromRef(ref,SamtoolsPath);
        Iterator<Map.Entry<String,Variant>> iter = this.vars.entrySet().iterator();

        VCFWriter vw = new VCFWriter(new File(this.outputFolder+"MaxEntScan_Filtered.vcf"),new File(this.ref+"hg19.fa"));
        VCFWriter sig_vw = new VCFWriter(new File(this.outputFolder+"MaxEntScan_Significant.vcf"),new File(this.ref+"hg19.fa"));
        VCFWriter possiblySig_vw = new VCFWriter(new File(this.outputFolder+"MaxEntScan_PossiblySignificant.vcf"),new File(this.ref+"hg19.fa"));
        VCFWriter notSig_vw = new VCFWriter(new File(this.outputFolder+"MaxEntScan_NonSignificant.vcf"),new File(this.ref+"hg19.fa"));

        while(iter.hasNext()){
            Map.Entry<String,Variant> entry = iter.next();
            Variant var = entry.getValue();
            var.parseSpliceInfo(rsp, prfr);
            if(this.algorithm.equals("MES")){
                MESRunner mr = new MESRunner(var,this.outputFolder,vw, this.algoritmPath);
                if (!mr.IsEmpty()){
                    int sig = var.checkMesSignificance();
                    if(sig == 2){
                        var.makeModifiedProtein();
                        System.out.println(Utilities.GREEN+"significant difference"+ Utilities.RESET);
                        VariantContextBuilder vcb = var.creatVariantContext();
                        vcb.attribute("MesScore",mr.getScores());
                        sig_vw.writeVar(vcb.make());
                    }
                    else if(sig == 1){
                        System.out.println(Utilities.GREEN+"possibly significant difference"+ Utilities.RESET);
                        VariantContextBuilder vcb = var.creatVariantContext();
                        vcb.attribute("MesScore",mr.getScores());
                        possiblySig_vw.writeVar(vcb.make());
                        iter.remove();
                        continue;
                    }
                    else if(sig == 0){
                        System.out.println(Utilities.GREEN+"not significant difference"+ Utilities.RESET);
                        VariantContextBuilder vcb = var.creatVariantContext();
                        vcb.attribute("MesScore",mr.getScores());
                        notSig_vw.writeVar(vcb.make());
                        iter.remove();
                        continue;
                    }
                    else{
                        iter.remove();
                        continue;
                    }

                }
            }


        }

        sig_vw.close();
        notSig_vw.close();
        possiblySig_vw.close();
        vw.close();
    }

    private void runAnnotations(String newFile) {
        AnnovarRunner AR = new AnnovarRunner(this.annovar, this.outputFolder);

        String gerp = AR.Gerp2(newFile, this.human);
        GeneralAnnotationParser parser = new GeneralAnnotationParser(this.outputFolder+gerp, true);
        this.vars = parser.parse(this.vars);

//        String phastCons = AR.PhastCons(newFile, this.human);
//        parser = new GeneralAnnotationParser(this.outputFolder+phastCons, false);
//        this.vars = parser.parsePhastCons(this.vars);

        String oneKGenomes = AR.onekGenomes(newFile,this.human);
        parser = new GeneralAnnotationParser(this.outputFolder+oneKGenomes,true);
        this.vars = parser.parse(this.vars);

        String exac = AR.Exac(newFile,this.human);
        parser = new GeneralAnnotationParser(this.outputFolder+exac, true);
        this.vars = parser.parse(this.vars);
    }

    private String convertAnnovar(){
        AnnovarRunner AR = new AnnovarRunner(this.annovar,this.outputFolder);
        String avinput = AR.convert2Annovar(this.input);
        String varfunc = AR.Gene(avinput,this.human);
        return varfunc;
    }

    private String writeNewAvinput() {
        try {
            String filename = "SpliceOnly.avinput";
            File newavinput = new File(this.outputFolder + filename);
            FileWriter writer = new FileWriter(newavinput);
            for (Map.Entry<String, Variant> entry : this.vars.entrySet()) {
                String key = entry.getKey();
                String[] chrpos = key.split(":");
                Variant v = entry.getValue();

                writer.write(chrpos[0] + "\t" + chrpos[1] + "\t" + chrpos[1] + "\t" + v.getRef() + "\t" + v.getAlt() + "\n");
            }
            writer.close();

            return filename;

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    private String geneBasedAnnotation(String varfunc){
        annovarWriter nonSplice = new annovarWriter(this.outputFolder+"NonSplice.txt");
        GeneParser gp = new GeneParser(this.outputFolder+varfunc);
        this.vars = gp.parse(nonSplice);
        return writeNewAvinput();
    }
}
