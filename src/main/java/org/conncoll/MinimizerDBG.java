package org.conncoll;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A "minimizer de Bruijn graph" builder, enhanced with advanced features for
 * error correction and graph cleaning, and adapted to work with the project's Menu class.
 *
 * This version retains the original method signatures called by the Menu, but
 * implements an advanced workflow internally.
 */
public class MinimizerDBG {

    // --- Parameters from Menu ---
    private final int k;
    private final int w;
    private final boolean verbose;

    // --- NEW: Hard-coded parameters for advanced features ---
    // These are set here because the Menu does not provide commands for them.
    private final int minAbundance = 2;
    private final int maxTipLength;

    // --- Core Data Structures ---
    private final Map<String, String> bucketedKmers = new HashMap<>();
    private final Map<String, Set<String>> adj = new HashMap<>();

    /**
     * Constructor called by the Menu class.
     */
    public MinimizerDBG(int k, int w, boolean verbose) {
        if (w > k) {
            throw new IllegalArgumentException("Window size w cannot exceed k-mer length k");
        }
        if (k <= 0 || w <= 0) {
            // Handle potentially uninitialized k/w values from Menu
            throw new IllegalArgumentException("k-mer and window size must be positive. Set them via the 'Kmer' and 'Window' commands.");
        }
        this.k = k;
        this.w = w;
        this.verbose = verbose;
        // Set the default tip length based on k
        this.maxTipLength = k * 2;
    }

    // --- UTILITY METHODS ---
    private String reverseComplement(String dna) {
        StringBuilder sb = new StringBuilder();
        for (int i = dna.length() - 1; i >= 0; i--) {
            char base = dna.charAt(i);
            switch (base) {
                case 'A': sb.append('T'); break;
                case 'T': sb.append('A'); break;
                case 'C': sb.append('G'); break;
                case 'G': sb.append('C'); break;
                default: sb.append('N'); break;
            }
        }
        return sb.toString();
    }

    private String getCanonical(String kmer) {
        String rcKmer = reverseComplement(kmer);
        return kmer.compareTo(rcKmer) <= 0 ? kmer : rcKmer;
    }

    private String computeMinimizer(String kmer) {
        String best = null;
        for (int i = 0; i <= k - w; i++) {
            String sub = kmer.substring(i, i + w);
            if (best == null || sub.compareTo(best) < 0) {
                best = sub;
            }
        }
        return best != null ? best.intern() : null;
    }

    // --- PUBLIC METHODS CALLED BY MENU ---

    /**
     * This is the main entry point called by the Menu. It now orchestrates the
     * entire advanced workflow: counting, filtering, building, and trimming.
     * @param fastqFile The input file provided from the Menu.
     */
    public void buildFromFastq(File fastqFile) throws IOException {
        if (verbose) {
            System.out.println("--- Starting Advanced DBG Construction Workflow ---");
            System.out.printf("Parameters: k=%d, w=%d, min_abundance=%d, max_tip_length=%d%n", k, w, minAbundance, maxTipLength);
        }

        // Pass 1: Count k-mer abundances
        Map<String, Integer> kmerCounts = countKmerAbundances(fastqFile);

        // Filter for "solid" k-mers
        Set<String> solidKmers = kmerCounts.entrySet().stream()
                .filter(entry -> entry.getValue() >= this.minAbundance)
                .map(Map.Entry::getKey)
                .collect(Collectors.toSet());
        kmerCounts = null; // Allow garbage collector to free memory
        if (verbose) System.out.printf("Identified %,d solid k-mers (abundance >= %d).%n", solidKmers.size(), this.minAbundance);

        // Pass 2: Build the graph using only solid k-mers
        buildFromSolidKmers(fastqFile, solidKmers);

        // Graph Cleaning Step
        trimTips();

        if (verbose) System.out.println("--- DBG Construction Workflow Complete ---");
    }

    /**
     * This method is called by the Menu after the graph is built. It operates on
     * the cleaned graph produced by the new `buildFromFastq` workflow.
     */
    public void writeUnitigsAsFasta(File outFasta) throws IOException {
        if (verbose) System.out.println("Writing unitigs to FASTA file: " + outFasta.getAbsolutePath());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFasta))) {
            long contigCount = 0;
            Set<String> visited = new HashSet<>();

            for (String minim : bucketedKmers.keySet()) {
                if (visited.contains(minim)) continue;
                int deg = adj.getOrDefault(minim, Collections.emptySet()).size();
                if (deg == 2) continue;

                List<String> path = walkFrom(minim, visited);
                String contigSeq = buildSequenceFromPath(path);
                if (contigSeq.length() >= k) {
                    contigCount++;
                    writer.write(String.format(">contig_%d_len_%d%n%s%n", contigCount, contigSeq.length(), contigSeq));
                }
            }

            for (String minim : bucketedKmers.keySet()) {
                if (visited.contains(minim)) continue;
                List<String> cyclePath = walkCycle(minim, visited);
                String contigSeq = buildSequenceFromPath(cyclePath);
                if (contigSeq.length() >= k) {
                    contigCount++;
                    writer.write(String.format(">contig_%d_cycle_len_%d%n%s%n", contigCount, contigSeq.length(), contigSeq));
                }
            }
            if (verbose) System.out.printf("Wrote %,d unitigs.%n", contigCount);
        }
    }


    // --- INTERNAL WORKFLOW METHODS (NOW PRIVATE) ---

    private Map<String, Integer> countKmerAbundances(File fastqFile) throws IOException {
        if (verbose) System.out.println("Pass 1: Counting k-mer abundances...");
        Map<String, Integer> kmerCounts = new HashMap<>();
        try (FastqReader reader = new FastqReader(fastqFile)) {
            for (FastqRecord rec : reader) {
                String seq = rec.getReadString().toUpperCase();
                if (seq.length() < k) continue;
                for (int i = 0; i <= seq.length() - k; i++) {
                    String kmer = seq.substring(i, i + k);
                    if (kmer.contains("N")) continue;
                    String canonicalKmer = getCanonical(kmer).intern();
                    kmerCounts.put(canonicalKmer, kmerCounts.getOrDefault(canonicalKmer, 0) + 1);
                }
            }
        }
        return kmerCounts;
    }

    private void buildFromSolidKmers(File fastqFile, Set<String> solidKmers) throws IOException {
        if (verbose) System.out.println("Pass 2: Building graph from solid k-mers...");
        try (FastqReader reader = new FastqReader(fastqFile)) {
            for (FastqRecord rec : reader) {
                String seq = rec.getReadString().toUpperCase();
                if (seq.length() < k) continue;
                String prevMin = null;
                for (int i = 0; i <= seq.length() - k; i++) {
                    String kmer = seq.substring(i, i + k);
                    if (kmer.contains("N")) {
                        prevMin = null;
                        continue;
                    }
                    String canonicalKmer = getCanonical(kmer).intern();
                    if (!solidKmers.contains(canonicalKmer)) {
                        prevMin = null;
                        continue;
                    }
                    String minim = computeMinimizer(canonicalKmer);
                    bucketedKmers.putIfAbsent(minim, canonicalKmer);
                    if (prevMin != null && !prevMin.equals(minim)) {
                        adj.computeIfAbsent(prevMin, x -> new HashSet<>()).add(minim);
                        adj.computeIfAbsent(minim, x -> new HashSet<>()).add(prevMin);
                    }
                    prevMin = minim;
                }
            }
        }
    }

    private void trimTips() {
        if (verbose) System.out.printf("Graph Cleaning: Trimming tips (paths <= %d nodes)...%n", this.maxTipLength);
        boolean changedInRound = true;
        int totalTipsRemoved = 0;

        while (changedInRound) {
            changedInRound = false;
            List<String> nodesToRemove = new ArrayList<>();
            for (String node : adj.keySet()) {
                if (adj.get(node).size() == 1) {
                    List<String> tipPath = new ArrayList<>();
                    String curr = node;
                    while (tipPath.size() <= this.maxTipLength) {
                        tipPath.add(curr);
                        Set<String> neighbors = adj.get(curr);
                        if (neighbors == null || neighbors.size() != 1) break;
                        curr = neighbors.iterator().next();
                    }
                    if (tipPath.size() <= this.maxTipLength) {
                        nodesToRemove.addAll(tipPath);
                    }
                }
            }
            if (!nodesToRemove.isEmpty()) {
                for (String node : new HashSet<>(nodesToRemove)) {
                    Set<String> neighbors = adj.remove(node);
                    if (neighbors != null) {
                        for (String neighbor : neighbors) {
                            Set<String> neighborAdj = adj.get(neighbor);
                            if (neighborAdj != null) neighborAdj.remove(node);
                        }
                        bucketedKmers.remove(node);
                        changedInRound = true;
                        totalTipsRemoved++;
                    }
                }
            }
        }
        if (verbose) System.out.printf("Finished trimming. Removed a total of %,d tip nodes.%n", totalTipsRemoved);
    }

    // --- PRIVATE GRAPH TRAVERSAL AND ASSEMBLY LOGIC ---

    private List<String> walkFrom(String startMinim, Set<String> visited) {
        List<String> path = new ArrayList<>();
        path.add(startMinim);
        visited.add(startMinim);

        // Walk “forward”
        String curr = startMinim;
        String prev = null;
        while (true) {
            Set<String> neighbors = adj.getOrDefault(curr, Collections.emptySet());
            if (adj.getOrDefault(curr, Collections.emptySet()).size() != 2 && !curr.equals(startMinim)) break;

            String next = null;
            for (String n : neighbors) {
                if (!n.equals(prev)) { next = n; break; }
            }
            if (next == null || visited.contains(next)) break;

            path.add(next); visited.add(next); prev = curr; curr = next;
        }

        // Walk “backward”
        List<String> backwardPath = new ArrayList<>();
        curr = startMinim;
        prev = path.size() > 1 ? path.get(1) : null;

        while (true) {
            Set<String> neighbors = adj.getOrDefault(curr, Collections.emptySet());
            if (adj.getOrDefault(curr, Collections.emptySet()).size() != 2 && !curr.equals(startMinim)) break;

            String next = null;
            for (String n : neighbors) {
                if (!n.equals(prev)) { next = n; break; }
            }
            if (next == null || visited.contains(next)) break;

            backwardPath.add(next); visited.add(next); prev = curr; curr = next;
        }

        Collections.reverse(backwardPath);
        backwardPath.addAll(path);
        return backwardPath;
    }

    private List<String> walkCycle(String minim, Set<String> visited) {
        List<String> cycle = new ArrayList<>();
        String start = minim;
        String curr = minim;
        String prev = null;
        do {
            cycle.add(curr);
            visited.add(curr);
            String next = null;
            for (String n : adj.getOrDefault(curr, Collections.emptySet())) {
                if (!n.equals(prev)) { next = n; break; }
            }
            prev = curr;
            curr = next;
        } while (curr != null && !curr.equals(start));
        return cycle;
    }

    private String buildSequenceFromPath(List<String> path) {
        if (path == null || path.isEmpty()) return "";
        String firstKmer = bucketedKmers.get(path.get(0));
        if (firstKmer == null) return "";

        StringBuilder contig = new StringBuilder(firstKmer);
        for (int i = 1; i < path.size(); i++) {
            String candidate = bucketedKmers.get(path.get(i));
            if (candidate == null) continue;
            contig.append(candidate.charAt(k - 1));
        }
        return contig.toString();
    }
}