package org.conncoll;

import com.google.common.base.Preconditions;
import com.google.common.hash.BloomFilter;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import org.conncoll.MDBG.Assembly;
import org.conncoll.MDBG.MDBGutils;
import org.conncoll.MDBG.MinimizerOccurrence;
import com.google.code.externalsorting.ExternalSort;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;


import java.io.*;


import static org.conncoll.MDBG.MDBGutils.*;

public class Menu {

    record Edge(int fromId, int toId, String sequence) {}
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
        outputFile = new File(path+"CONTIGS_"+inputFile.getName());
    }

    public static void setSam(String path){
        samFile = new File(path);
    }

    public static void setFilteredFastqOut(String path){

        filteredFastqFileOut = new File(path+"FILTERED_"+inputFile.getName());
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
    public static double density = .2;


    public static Long2IntOpenHashMap minimizerToId = new Long2IntOpenHashMap();
    public static int nextMinimizerId = 0;


    public static void build() throws IOException {
        // --- SETUP ---
        Preconditions.checkArgument(kmerSize > 0, "Kmer size must be larger than 0");
        Preconditions.checkArgument(inputFile.exists(), "Invalid input file path");
        String unsortedEdgePath = "/Volumes/T9/Raw_Data/JMDBGout/EdgeStoreage/Edges_Unsorted.tmp";
        String sortedEdgePath = "/Volumes/T9/Raw_Data/JMDBGout/EdgeStoreage/Edges_Sorted.tmp";
        File fastqToProcess;

        // Re-initialize for each build
        minimizerToId = new Long2IntOpenHashMap();
        nextMinimizerId = 0;

        // ---  STEP 1: PRE-FILTERING (Optional) ---
        if (samFile != null && samFile.exists()) {
            System.out.println("Pre-filtering reads using SAM file...");
            Set<String> matchedReads = samUtilities.getMappedReadNames(samFile);
            fastqToProcess = FastQFilterer.filterFastq(inputFile, filteredFastqFileOut, matchedReads);
        } else {
            System.out.println("SKIPPING SAM PRE-FILTERING: No SAM file provided.");
            fastqToProcess = inputFile;
        }

        // ---  STEP 2: BLOOM FILTERS (RESTORED and ENABLED) ---
        System.out.println("Building Bloom filters to find frequent k-mers...");
        BloomFilter<Long> seenOnce = multiThreadedBF(fastqToProcess, 7, null); // Using 7 threads as an example
        BloomFilter<Long> seenTwice = multiThreadedBF(fastqToProcess, 7, seenOnce);
        System.out.println("Bloom filters built.");


        // ---  STEP 3 (PASS 1): DISCOVER ALL MINIMIZER NODES ---
        System.out.println("PASS 1: Discovering all minimizer nodes...");
        long threshold = (long) (density * (double) Long.MAX_VALUE);
        try (BufferedReader reader = new BufferedReader(new FileReader(fastqToProcess))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = reader.readLine(); // Sequence
                reader.readLine(); // +
                reader.readLine(); // Quality
                if (line == null) break;

                EncodedSequence encoded = encode(line);
                int validBaseCount = encoded.validBaseCount();
                if (validBaseCount < kmerSize) continue;

                byte[] encodedSequence = encoded.sequence();
                long KMER_MASK = (1L << (kmerSize * 2)) - 1;
                long currentKmer = 0L;

                for (int i = 0; i < validBaseCount; i++) {
                    currentKmer <<= 2;
                    currentKmer |= getBaseAt(encodedSequence, i);
                    if (i >= kmerSize - 1) {
                        currentKmer &= KMER_MASK;
                        long canonicalKmer = MDBGutils.getCanonical(currentKmer, kmerSize);

                        if (seenTwice.mightContain(canonicalKmer)) {
                            long h = fnv1a64(canonicalKmer);
                            if ((h & Long.MAX_VALUE) < threshold) {
                                minimizerToId.computeIfAbsent(canonicalKmer, key -> nextMinimizerId++);
                            }
                        }
                    }
                }
            }
        }
        System.out.println("PASS 1 Complete. Found " + minimizerToId.size() + " unique minimizer nodes.");

        // ---  STEP 4 (PASS 2): BUILD EDGES FROM READS ---
        System.out.println("PASS 2: Building edges from reads...");
        Set<Edge> uniqueEdges = new HashSet<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(fastqToProcess))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = reader.readLine(); // Sequence
                reader.readLine(); // + stuff
                reader.readLine(); // Quality
                if (line == null) break;

                List<MinimizerOccurrence> occurrencesInRead = new ArrayList<>();
                EncodedSequence encoded = encode(line);
                int validBaseCount = encoded.validBaseCount();
                if (validBaseCount < kmerSize) continue;

                byte[] encodedSequence = encoded.sequence();
                long KMER_MASK = (1L << (kmerSize * 2)) - 1;
                long currentKmer = 0L;

                for (int i = 0; i < validBaseCount; i++) {
                    currentKmer <<= 2;
                    currentKmer |= getBaseAt(encodedSequence, i);
                    if (i >= kmerSize - 1) {
                        currentKmer &= KMER_MASK;
                        long canonicalKmer = MDBGutils.getCanonical(currentKmer, kmerSize);
                        if (minimizerToId.containsKey(canonicalKmer)) {
                            int id = minimizerToId.get(canonicalKmer);
                            int pos = i - kmerSize + 1;
                            occurrencesInRead.add(new MinimizerOccurrence(id, pos));
                        }
                    }
                }

                for (int i = 0; i < occurrencesInRead.size() - 1; i++) {
                    MinimizerOccurrence fromMinimizer = occurrencesInRead.get(i);
                    MinimizerOccurrence toMinimizer = occurrencesInRead.get(i + 1);
                    if(fromMinimizer.id == toMinimizer.id) continue;
                    String sequence = line.substring(fromMinimizer.pos, toMinimizer.pos + kmerSize);
                    uniqueEdges.add(new Edge(fromMinimizer.id, toMinimizer.id, sequence));
                }
            }
        }

        Int2IntOpenHashMap inDegrees = new Int2IntOpenHashMap();
        Int2IntOpenHashMap outDegrees = new Int2IntOpenHashMap();
        for (Edge edge : uniqueEdges) {
            outDegrees.addTo(edge.fromId(), 1);
            inDegrees.addTo(edge.toId(), 1);
        }
        System.out.println("PASS 2 Complete. Generated " + uniqueEdges.size() + " unique edges.");

        // --- SORTING AND ASSEMBLY ---
        System.out.println("Writing and sorting " + uniqueEdges.size() + " edges...");
        try (BufferedWriter edgeWriter = new BufferedWriter(new FileWriter(unsortedEdgePath))) {
            for (Edge edge : uniqueEdges) {
                edgeWriter.write(edge.fromId() + "\t" + edge.toId() + "\t" + edge.sequence() + "\n");
            }
        }

        Comparator<String> comparator = (line1, line2) -> {
            String fromId1 = line1.substring(0, line1.indexOf('\t'));
            String fromId2 = line2.substring(0, line2.indexOf('\t'));
            return Integer.compare(Integer.parseInt(fromId1), Integer.parseInt(fromId2));
        };

        try {
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
