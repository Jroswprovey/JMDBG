package org.conncoll;

import com.google.common.base.Preconditions;
import com.google.common.hash.BloomFilter;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.conncoll.MDBG.Assembly;
import org.conncoll.MDBG.MinimizerOccurrence;
import com.google.code.externalsorting.ExternalSort;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Comparator;


import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;


import static org.conncoll.MDBG.MDBGutils.*;

public class Menu {

    public static File inputFile;
    public static File outputFile;
    public static File samFile;
    public static File filteredFastqFileOut;
    public static boolean verbose = false;
    public static boolean compress;
    public static int kmerSize;


    public static List<Command> commands = List.of(

            new Command("Help", "Lists commands", Menu::helpList),
            new Command("InputF", "should be followed with an input fastq file", Menu::handleInputFile),
            new Command("OutputF", "Define a out path (must end in .fastq)", Menu::handleOutputFile),
            new Command("SamF", "Define a sam file for filtering new fastq file", Menu::setSam),
            new Command("SamOut", "Set the output location for the filtered fastq file", Menu::setFilteredFastqOut),
            new Command("Verbose", "Outputs debug info while running", Menu::setVerbose),
            new Command("FQRead", "Reads a fastq file, ", Menu::readFastQ),
            new Command("Compress","Followed by true of false to idicate weather the output file will be compressed (Defaults to the input files state)",  Menu::setCompress),
            new Command("Kmer",    "set k-mer length (integer)",  Menu::setKmerSize),
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
        outputFile = new File(path+"CONTIGS"+inputFile.getName());
    }

    public static void setSam(String path){
        samFile = new File(path);
    }

    public static void setFilteredFastqOut(String path){

        filteredFastqFileOut = new File(path+"FILTERED"+inputFile.getName());
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

    public static void setKmerSize(String arg){
        kmerSize = Integer.parseInt(arg);
    }

    public static volatile long  totalKmers = 0;
    public static volatile long totalMinimizers = 0;
    public static double density = .005;


    public static Long2IntOpenHashMap minimizerToId = new Long2IntOpenHashMap();
    public static int nextMinimizerId = 0;

    public static void build() throws IOException {
        //checking users inputs
        Preconditions.checkArgument(kmerSize > 0, "Kmer size must be larger than 0");
        Preconditions.checkArgument(inputFile.exists(), "Invalid input file path");
        Preconditions.checkArgument(samFile.exists(), "Invalid sam file path");



        //edges are stored in a temporary location on disk to ease pressure on memory
        String unsortedEdgePath = "/Volumes/T9/Raw_Data/JMDBGout/EdgeStoreage/Edges_Unsorted.tmp";
        String sortedEdgePath = "/Volumes/T9/Raw_Data/JMDBGout/EdgeStoreage/Edges_Sorted.tmp";

        Int2IntOpenHashMap inDegrees = new Int2IntOpenHashMap();
        Int2IntOpenHashMap outDegrees = new Int2IntOpenHashMap();

        //this is run on a separate thread to not bog down the actual work with prints to console
        ProgressDisplay progressDisplay = new ProgressDisplay();
        Thread displayThread = new Thread(progressDisplay);
        displayThread.start();

        Set<String> matchedReads = samUtilities.getMappedReadNames(samFile);


        //Step 1: Pre filter file with aligner
        File filteredFastQ = FastQFilterer.filterFastq(inputFile, filteredFastqFileOut,matchedReads);


        //Step 2: filter file with bloom filters
        BloomFilter seenOnce = multiThreadedBF(filteredFastQ, 7, null);
        BloomFilter seenTwice = multiThreadedBF(filteredFastQ, 7, seenOnce);

        saveBloomFilter(seenTwice, "/Volumes/T9/Raw_Data/JMDBGout/BloomFilters/seenTwice.bf");
        //BloomFilter seenTwice = readBloomFilter("/Volumes/T9/Raw_Data/JMDBGout/BloomFilters/seenTwice.bf");



        //Step 3: select kmers based on a threshold

        try(BufferedReader reader = new BufferedReader(new FileReader(filteredFastQ));
            BufferedWriter edgeWriter = new BufferedWriter(new FileWriter(unsortedEdgePath))){
            String line;
            long threshold = (long)(density * (double)Long.MAX_VALUE);
                while (reader.readLine() !=null) {//ignore header

                    line = reader.readLine(); // actual sequence only store this
                    reader.readLine(); // + stuff
                    reader.readLine(); // quality scores

                    List<MinimizerOccurrence> occurrencesInRead = new ArrayList<>(); //resets every sequence so memory footprint stays low



                    EncodedSequence encoded = encode(line); //encode the sequence read in from earlier
                    int validBaseCount = encoded.validBaseCount();
                    byte[] encodedSequence = encoded.sequence();


                    if (validBaseCount < kmerSize) {
                        continue; // Sequence is shorter than k.
                    }

                    // A mask to keep our k-mer at 62 bits (31 bases * 2 bits/base).
                    long KMER_MASK = (1L << (kmerSize * 2)) - 1;

                    // --- Step 1: Manually build the first k-mer ---
                    long currentKmer = 0L;
                    for (int j = 0; j < kmerSize; j++) {
                        currentKmer <<= 2; // Make room for the next base.
                        byte newBaseCode = getBaseAt(encodedSequence, j); // Get 2-bit code for base i.
                        currentKmer |= newBaseCode; // Add it to the k-mer.
                    }

                    if (seenTwice.mightContain(currentKmer)) {

                        Menu.totalKmers++;

                        //hash the current kmer
                        long h = fnv1a64(currentKmer);
                        long unsigned = h & Long.MAX_VALUE;// make non-negative


                        //if the current kmer falls below the threshold, add
                        if (unsigned < threshold) {
                            Menu.totalMinimizers++;
                            //if the string lookup doesn't contain it
                                int id = minimizerToId.computeIfAbsent(currentKmer, key -> nextMinimizerId++);
                                occurrencesInRead.add(new MinimizerOccurrence(id, 0));


                        }
                    }


                    // --- Step 2: Slide the window for the rest of the sequence ---
                    for (int k = kmerSize; k < validBaseCount; k++) {
                        // a. Shift the k-mer left by 2 bits, discarding the oldest base.
                        currentKmer <<= 2;

                        // b. Get the 2-bit code for the new base entering the window.
                        byte newBaseCode = getBaseAt(encodedSequence, k);

                        // c. Add the new base to the right side of the k-mer.
                        currentKmer |= newBaseCode;

                        // d. Apply the mask to trim the k-mer back to the correct 62-bit length.
                        currentKmer &= KMER_MASK;

                        // Process the new overlapping k-mer
                        if (seenTwice.mightContain(currentKmer)) {
                            Menu.totalKmers++;

                            //hash the current kmer
                            long h = fnv1a64(currentKmer);
                            long unsigned = h & Long.MAX_VALUE;// make non-negative

                            //if the current kmer falls below the threshold, add
                            if (unsigned < threshold) {
                                Menu.totalMinimizers++;
                                //if the string lookup doesn't contain it
                                    int id = minimizerToId.computeIfAbsent(currentKmer, key -> nextMinimizerId++);
                                    occurrencesInRead.add(new MinimizerOccurrence(id, k - kmerSize + 1));

                            }
                        }
                    }

                    for (int i = 0; i < occurrencesInRead.size() - 1; i++) {
                        MinimizerOccurrence fromMinimizer = occurrencesInRead.get(i);
                        MinimizerOccurrence toMinimizer = occurrencesInRead.get(i + 1);

                        outDegrees.addTo(fromMinimizer.id, 1);
                        inDegrees.addTo(toMinimizer.id, 1);

                        int start = fromMinimizer.pos;
                        int end = toMinimizer.pos + kmerSize;
                        if (end > line.length()) {
                            end = line.length();
                        }

                        String sequence = line.substring(start, end);
                        //write minimizers to disk (may take very long :/ possible multithreading in the future,
                        //but current design isn't friendly to that)
                        edgeWriter.write(
                                fromMinimizer.id + "\t" +
                                        toMinimizer.id + "\t" +
                                        sequence + "\n"
                        );


                    }
                }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.gc();

        System.out.println("Sorting edges using Java external sort...");

        Comparator<String> comparator = (line1, line2) -> {
            String fromId1 = line1.substring(0, line1.indexOf('\t'));
            String fromId2 = line2.substring(0, line2.indexOf('\t'));
            return Integer.compare(Integer.parseInt(fromId1), Integer.parseInt(fromId2));
        };

        try {
            // This one line performs the entire external sort
            ExternalSort.mergeSortedFiles(
                    ExternalSort.sortInBatch(new File(unsortedEdgePath), comparator),
                    new File(sortedEdgePath),
                    comparator,
                    StandardCharsets.UTF_8
            );
            System.out.println("Java external sort completed successfully.");
        } catch (IOException e) {
            throw new RuntimeException("Error during Java external sort", e);
        }

        new Assembly().assembleStreaming(sortedEdgePath, outputFile.getAbsolutePath(), kmerSize, inDegrees, outDegrees);

        progressDisplay.isFinished = true;
        try {
            displayThread.join();
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
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
