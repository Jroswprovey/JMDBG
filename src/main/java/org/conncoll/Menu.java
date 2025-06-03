package org.conncoll;

import java.io.File;
import java.util.List;

public class Menu {

    public static File inputFile;
    public static File outputFile;
    public static File samFile;


    public static List<Command> commands = List.of(

            new Command("Help", "Lists commands", Menu::helpList),
            new Command("InputF", "should be followed with an input fastq file", (arg) -> {handleInputFile(arg);}),
            new Command("OutputF", "Define a out path (must end in .fastq)", (arg) -> {handleOutputFile(arg);}),
            new Command("SamF", "Define a sam file for filtering new fastq file", (arg) -> {setSam(arg);})
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






}
