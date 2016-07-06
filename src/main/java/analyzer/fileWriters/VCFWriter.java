package analyzer.fileWriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.EnumSet;
import java.util.HashSet;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

/**
 * Created by mwads on 1/27/16.
 */
public class VCFWriter {


    private VariantContextWriter writer;
    private String path;

    public VCFWriter(File file, File RefPath){
        try {
            SAMSequenceDictionary dict = new IndexedFastaSequenceFile(RefPath).getSequenceDictionary();
            VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
            HashSet<VCFHeaderLine> headerlines = new HashSet<>();
            headerlines.add(new VCFHeaderLine("INFO","<ID=Gene,Number=1,Type=String,Description=The Gene name>"));
            headerlines.add(new VCFHeaderLine("INFO","<ID=Transcripts,Number=1,Type=String,Description=The transcript names of variants that were more than 20 bases off the 3' end of the intron and 6 bases off the 5' end of the intron>"));
            headerlines.add(new VCFHeaderLine("INFO","<ID=MesScore,Number=1,Type=Double,Description=The MaxEntScan score for the variant>"));

            VCFHeader vh = new VCFHeader(headerlines);

            builder = builder.setOutputFile(file);

            builder = builder.setReferenceDictionary(dict);

            EnumSet<Options> ops = EnumSet.of(Options.INDEX_ON_THE_FLY,Options.ALLOW_MISSING_FIELDS_IN_HEADER,Options.DO_NOT_WRITE_GENOTYPES);
            builder = builder.setOptions(ops);

            this.writer = builder.build();
            this.writer.writeHeader(vh);
            this.path = file.getAbsolutePath();



        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }

    public void writeVar(VariantContext vc){
        System.out.println(vc.toString());
        this.writer.add(vc);
    }

    public String getPath(){
        return this.path;
    }

    public void close() {
        this.writer.close();
    }
}
