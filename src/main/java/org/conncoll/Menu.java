package org.conncoll;

import java.io.File;
import java.util.List;

public class Menu {



    public static List<Command> commands = List.of(

            new Command("Help", "Lists commands", Menu::helpList),
            new Command("InputF", "should be followed with an input fastq file", (arg) -> {handleInputFile(arg);} ),
            new Command("OutputF", "Define a out path (must end in .fastq)", Menu::handleOutputFile)


    );

    public static void helpList(){
        for (Command command : commands) {
            System.out.println(command.getName() + ": " + command.getDescription());

        }
    }

    public static File handleInputFile(String path){
        return new File(path);
    }

    public static File handleOutputFile(String path){
        return new File(path);
    }


}
