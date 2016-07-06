package analyzer.AnnovarRunners;

import analyzer.Utilities.Utilities;
import net.sourceforge.argparse4j.inf.Namespace;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;

/**
 * Created by mwads on 1/15/16.
 */
public class AnnovarRunner {

    private String AnnovarPath;
    private String OutputFolder;





    public AnnovarRunner(String Path2Annovar, String OutputFolder){
        this.AnnovarPath = Path2Annovar;
        this.OutputFolder = OutputFolder;
    }

    public String convert2Annovar(String vcf) {
        System.out.println(Utilities.GREEN+"Running VCF to Annovar Conversion"+ Utilities.RESET);
        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/convert2annovar.pl","-format","vcf4",vcf,"--outfile",vcf+".avinput"};
            ProcessBuilder pb = new ProcessBuilder(call);

            Process p = pb.start();
            p.waitFor();
            System.out.println(Utilities.getProcessOutput(p));
            System.out.println(Utilities.getProcessError(p));
            p.destroy();

            String[] temp = vcf.split("/");
            return temp[temp.length-1]+".avinput";
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;
    }

    public String Gene(String avinput, String human){
        System.out.println(Utilities.GREEN+"Running Gene Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","--splicing_threshold","50","-build","hg19","-hgvs","-out",this.OutputFolder+avinput,avinput,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            return avinput+".variant_function";

        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;

    }

    public String PhastCons(String input, String human){
        System.out.println(Utilities.GREEN+"Running PhastCons Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-regionanno","-build","hg19","-dbtype","phastConsElements46way","-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            return input+".hg19_phastConsElements46way";

        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;
    }

    public String onekGenomes(String input, String human){
        System.out.println(Utilities.GREEN+"Running 1000 Genomes MAF Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-filter","-dbtype","1000g2015aug_all","-build","hg19","-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            return input+".hg19_ALL.sites.2015_08_";
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;
    }

    public String Gerp2(String input, String human){
        System.out.println(Utilities.GREEN+"Running Gerp++ Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-filter","-dbtype","gerp++gt2","-build","hg19","-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            return input+".hg19_gerp++gt2_";
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;

    }

    public String Exac(String input, String human){
        System.out.println(Utilities.GREEN+"Running Exac Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-filter","-dbtype","exac03","-build","hg19","-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            return input+".hg19_exac03_";
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;

    }


    //  Runtime.getRuntime().exec()
    /*
         convert2annovar.pl -format vcf4 example/ex2.vcf > ex2.avinput
         annotate_variation.pl -out ex1 -build hg19 -hgvs example/ex1.avinput humandb/
        annotate_variation.pl -regionanno -build hg19 -out ex1 -dbtype phastConsElements46way example/ex1.avinput humandb/ -normscore_threshold 400
        annotate_variation.pl -filter -dbtype 1000g2012apr_eur -buildver hg19 -out ex1 example/ex1.avinput humandb/
        annotate_variation.pl -filter -build hg19 -dbtype exac03 example/ex1.avinput humandb/
        annotate_variation.pl -filter -dbtype gerp++gt2 -out ex1 -build hg19 example/ex1.avinput humandb/

     */

}
