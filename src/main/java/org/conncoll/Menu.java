package org.conncoll;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.conncoll.MDBG.MinimizerDeBruijnGraph;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.List;
import java.util.concurrent.TimeUnit;

public class Menu {

    public static File inputFile;
    public static File outputFile;
    public static File samFile;
    public static boolean verbose = false;
    public static boolean compress;
    public static int kmer;
    public static int window;
    public static File minimizerInputFile;
    public static File minimizerOutputFile;

    public static List<Command> commands = List.of(

            new Command("Help", "Lists commands", Menu::helpList),
            new Command("InputF", "should be followed with an input fastq file", Menu::handleInputFile),
            new Command("OutputF", "Define a out path (must end in .fastq)", Menu::handleOutputFile),
            new Command("SamF", "Define a sam file for filtering new fastq file", Menu::setSam),
            new Command("Verbose", "Outputs debug info while running", Menu::setVerbose),
            new Command("FQRead", "Reads a fastq file, ", Menu::readFastQ),
            new Command("Compress","Followed by true of false to idicate weather the output file will be compressed (Defaults to the input files state)",  Menu::setCompress),
            new Command("Kmer",    "set k-mer length (integer)",  Menu::setKmer),
            new Command("Window",  "set window size for minimizers (integer)", Menu::setWindow),
            new Command("MinInput", "Path to FASTQ input file for minimizer graph", (arg) -> {Menu.minimizerInputFile = new File(arg);}),
            new Command("MinOutput", "Path to write unitigs FASTA output", (arg) -> {Menu.minimizerOutputFile = new File(arg);}),
            new Command("Filter","Filters a Fastq file given a sam file", Menu::filter),
            new Command("Assemble", "assemble MDBG", () -> {
                try {
                    build();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            })
    );

    public static void helpList(){
        for (Command command : commands) {
            System.out.println(command.getName() + ": " + command.getDescription());

        }
    }

    public static void handleInputFile(String path){
        inputFile = new File(path);
    }

    public static void handleOutputFile(String path){
        outputFile = new File(path);
    }

    public static void setSam(String path){
        samFile = new File(path);
    }

    public static void setVerbose(){
        verbose = true;
    }

    public static void readFastQ(String path){
        FQreader.setFile(new File(path));
        try {
            FQreader.readRecord();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    public static void setCompress(String arg){
        compress = Boolean.parseBoolean(arg);
    }

    public static void setKmer(String arg){
        kmer = Integer.parseInt(arg);
    }

    public static void setWindow(String arg){
        window = Integer.parseInt(arg);
    }

    public static void build() throws IOException {


    }

    public static void filter(){
        System.out.println(inputFile.getAbsolutePath());
        try {
            FastQFilterer.filterFastq(inputFile,outputFile,samUtilities.getMappedReadNames(samFile));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


}
