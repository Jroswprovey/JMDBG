package org.conncoll;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class Menu {

    public static File inputFile;
    public static File outputFile;
    public static File samFile;
    public static Boolean verbose = false;
    public static boolean compress;

    public static List<Command> commands = List.of(

            new Command("Help", "Lists commands", Menu::helpList),
            new Command("InputF", "should be followed with an input fastq file", Menu::handleInputFile),
            new Command("OutputF", "Define a out path (must end in .fastq)", Menu::handleOutputFile),
            new Command("SamF", "Define a sam file for filtering new fastq file", Menu::setSam),
            new Command("Verbose", "Outputs debug info while running", Menu::setVerbose),
            new Command("FQRead", "Reads a fastq file, ", Menu::readFastQ),
            new Command("Compess","Followed by true of false to idicate weather the output file will be compressed (Defaults to the input files state)",  Menu::setCompress)
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
            FQreader.read();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    public static void setCompress(String arg){
        compress = Boolean.parseBoolean(arg);
    }


}
