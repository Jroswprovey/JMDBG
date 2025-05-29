/*
Joseph Rosw-Provey
May 2025
This class sets out to handle communication with the minimap2 program.

 */

package org.conncoll;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

public class MinimapController {

    public static void analyze(String ref, String query, String destinationPath) throws IOException, InterruptedException {

        File outputFile = new File(destinationPath);
        ProcessBuilder pb = new ProcessBuilder("minimap2",ref,query);

        //settings for pb
        pb.redirectErrorStream(true);
        pb.redirectOutput(outputFile);

        Process process = pb.start();
        System.out.println("Starting MiniMap2 With command: " + "\n" + pb.command());


        try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                System.out.println(line); // Print any messages (errors, verbose output) to the console
                System.out.flush();       // Ensure it's printed immediately
            }
        }

       int exitcode = process.waitFor();
       System.out.println(exitcode);

    }


}
