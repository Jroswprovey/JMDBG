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

        long beginTime = System.nanoTime();


        System.out.println("Awesome Joe contig assembler");
        System.out.println("Setting up parameters...");
        int k = 21;
        int w = 31;

        File fastqFile = new File("/Volumes/T9/Raw Data/Sequence Read Files/PB644_EB813.hifi_reads.fastq.gz");
        File edgeStorageFile = new File("edge_sequences.tmp");

        if (!fastqFile.exists()) {
            System.err.println("ERROR: File not found at path: " + fastqFile.getAbsolutePath());
            return;
        }


        System.out.println("Pass 1: Counting total records in file...");
        long totalRecords = 0;
        try (FastqReader counterReader = new FastqReader(fastqFile)) {
            while (counterReader.hasNext()) {
                counterReader.next();
                totalRecords++;
            }
        }
        System.out.println("Found " + totalRecords + " total records. Starting main processing.");

        MinimizerDeBruijnGraph mdbg = new MinimizerDeBruijnGraph(k, w);
        long recordsProcessed = 0;
        long lastUpdateTime = System.nanoTime();
        long thirtySecondsInNanos = TimeUnit.SECONDS.toNanos(5);

        try (RandomAccessFile edgeFile = new RandomAccessFile(edgeStorageFile, "rw")) {
            try (FastqReader reader = new FastqReader(fastqFile)) {
                while (reader.hasNext()) {
                    FastqRecord record = reader.next();
                    mdbg.addSequence(record.getReadString(), edgeFile);
                    recordsProcessed++;

                    long currentTime = System.nanoTime();
                    if (currentTime - lastUpdateTime > thirtySecondsInNanos) {
                        double percentage = ((double) recordsProcessed / totalRecords) * 100.0;
                        System.out.printf("\rProgress: %.2f%% (%d / %d records processed)", percentage, recordsProcessed, totalRecords);
                        lastUpdateTime = currentTime;
                    }
                }
            }
        }

        System.out.printf("\rProgress: 100.00%% (%d / %d records processed)\n", totalRecords, totalRecords);
        System.out.println("Graph construction complete!");

        // --- 3. ASSEMBLE CONTIGS (With Corrections) ---
        System.out.println("\nAssembling contigs...");


        List<List<Long>> allContigPaths = mdbg.assembleContigs();
        System.out.println("Found " + allContigPaths.size() + " contigs.");

        int contigNumber = 1;

        for (List<Long> path : allContigPaths) {

            String dnaSequence = mdbg.stitchContig(path, edgeStorageFile);

            System.out.println("\n>Contig_" + contigNumber + " length=" + dnaSequence.length() + "bp");
            System.out.println(dnaSequence);
            contigNumber++;
        }


        edgeStorageFile.delete();

        long endTime = System.nanoTime();
        System.out.println("\nTotal execution time: " + TimeUnit.NANOSECONDS.toMinutes(endTime - beginTime) + " minutes.");
    }
}