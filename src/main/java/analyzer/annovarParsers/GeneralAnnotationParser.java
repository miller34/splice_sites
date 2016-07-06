package analyzer.annovarParsers;

import analyzer.variantInfo.Variant;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.TreeMap;

/**
 * Created by mwads on 1/15/16.
 */
public class GeneralAnnotationParser {

    private Scanner dropped;
    private Scanner filtered;
    private Scanner phastCons;


    public GeneralAnnotationParser(String annovarFile, Boolean twoFiles){
        try {

            if(twoFiles) {
                this.dropped = new Scanner(new File(annovarFile + "dropped"));
                this.filtered = new Scanner(new File(annovarFile + "filtered"));
            }
            else{
                this.phastCons = new Scanner(new File(annovarFile));
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public TreeMap<String, Variant> parse(TreeMap<String, Variant> Vars){
        while (this.dropped.hasNextLine()){
            String line = this.dropped.nextLine();
            String[] columns =line.split("\t");

            StringBuilder sb = new StringBuilder();
            sb.append(columns[2]+":"+columns[3]);

            String var = sb.toString();
            String anno = columns[1];
            Variant v = Vars.get(var);
            v.addAnnotation(anno);
            Vars.put(var, v);
        }
        while (this.filtered.hasNextLine()){
            String line = this.filtered.nextLine();
            String[] columns =line.split("\t");

            StringBuilder sb = new StringBuilder();
            sb.append(columns[0]+":"+columns[1]);

            String var = sb.toString();
            String anno = "NA";
            Variant v = Vars.get(var);
            v.addAnnotation(anno);
            Vars.put(var, v);
        }
        return Vars;
    }

     public TreeMap<String, Variant> parsePhastCons(TreeMap<String, Variant> Vars){
        while (this.phastCons.hasNextLine()){
            String line = this.phastCons.nextLine();
            String[] columns =line.split("\t");

            StringBuilder sb = new StringBuilder();
            sb.append(columns[2]+":"+columns[3]);

            String var = sb.toString();
            String anno = columns[1];
            Variant v = Vars.get(var);
            v.addAnnotation(anno);
            Vars.put(var, v);
        }
        return Vars;
    }

}
