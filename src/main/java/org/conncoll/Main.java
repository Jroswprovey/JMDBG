package org.conncoll;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.conncoll.MDBG.MinimizerDeBruijnGraph;

import java.io.RandomAccessFile;
import java.util.concurrent.TimeUnit;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Main {
    public static void main(String[] args) throws IOException {

//        File samFile = new File("/Volumes/T9/Raw Data/Aligned Sequence Files/PB644_EB816.hifi_reads.sam");
//        File bamFile = new File("/Volumes/T9/Raw Data/Bams/PB644_EB816.hifi_reads.bam");

        //for the length of the arg array
        for(int i =0; i < args.length; i++){

            //compare it to everything in the Menu list
            for(int j = 0; j < Menu.commands.size(); j++) {

                //check if a match is found
                if(args[i].equalsIgnoreCase(Menu.commands.get(j).getName())){
                    /*checks if the command being executed in a consumer, if it is, it requires a string to be passed
                     *along to be executed properly
                     */
                    if(Menu.commands.get(j).isConsumer()){
                        Menu.commands.get(j).execute(args[i+1]);
                    }else {
                        Menu.commands.get(j).execute();
                    }
                }
            }
        }
    }
}