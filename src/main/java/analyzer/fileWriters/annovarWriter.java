package analyzer.fileWriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

/**
 * Created by mwads on 1/29/16.
 */
public class annovarWriter {

    private FileWriter writer;

    public annovarWriter(String fileName){
        try {
            this.writer = new FileWriter(new File(fileName));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeLine(String line){
        try {
            this.writer.write(line+"\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void closeWriter(){
        try {
            this.writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
