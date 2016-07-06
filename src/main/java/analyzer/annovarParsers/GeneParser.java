package analyzer.annovarParsers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import analyzer.fileWriters.annovarWriter;
import analyzer.variantInfo.Variant;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.TreeMap;

/**
 * Created by mwads on 1/15/16.
 */
public class GeneParser {

    private Scanner varFunct;

    public GeneParser(String annovarGeneFile){

        try {
            this.varFunct = new Scanner(new File(annovarGeneFile));


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public TreeMap<String,Variant> parse(annovarWriter nonSplice){

        TreeMap<String, Variant> vars = new TreeMap<>();

        while(varFunct.hasNextLine()){
            String Line = varFunct.nextLine();
            String[] columns = Line.split("\t");
            if (columns[0].equals("splicing")){

                if (!columns[5].equals("-") && !columns[6].equals("-")) { //Excludes indels
                    Variant var = new Variant(columns[2], Integer.valueOf(columns[3]), columns[1], columns[5], columns[6], columns[7]);
                    StringBuilder key = new StringBuilder(columns[2]+":"+columns[3]);
                    vars.put(key.toString(), var);
                }
                else{ //indels get written to same file as others.
//                    VariantContext vc = buildContext(columns);
                    StringBuilder sb = new StringBuilder(columns[2]+"\t"+columns[3]+"\t"+columns[4]+"\t"+columns[5]+"\t"+columns[6]+"\t"+columns[7]+"\t"+columns[8]+"\t"+columns[9]+"\t"+columns[0]+"\t"+columns[1]);
                    nonSplice.writeLine(sb.toString());
                }
            }
            else{//These go to the vcf writer to write the non splice site file.
//                VariantContext vc = buildContext(columns);
                StringBuilder sb = new StringBuilder(columns[2]+"\t"+columns[3]+"\t"+columns[4]+"\t"+columns[5]+"\t"+columns[6]+"\t"+columns[7]+"\t"+columns[8]+"\t"+columns[9]+"\t"+columns[0]+"\t"+columns[1]);
                nonSplice.writeLine(sb.toString());
            }
        }
        return vars;
    }

    private VariantContext buildContext(String[] columns){
        System.out.println("In GeneParser");
        ArrayList<Allele> alleles = new ArrayList<>();
        if (!columns[5].equals("-"))
            alleles.add(Allele.create(columns[5],true));
        else
            alleles.add(Allele.create(".", true));

        if (!columns[6].equals("-")) {
            System.out.println(columns[6]);
            alleles.add(Allele.create(columns[6], false));
        }
        else
            alleles.add(Allele.create(".",false));

        VariantContextBuilder builder = new VariantContextBuilder("Annovar Exonic", columns[2], Long.valueOf(columns[3]), Long.valueOf(columns[4]), alleles);

        return builder.make();
    }
}

