package org.conncoll;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.*;
import java.util.*;

/**
 * A minimizer de Bruijn graph builder:
 * - Operates on canonical k-mers (lexicographically smallest of k-mer and its reverse complement).
 * - Buckets canonical k-mers by their w-length minimizer (computed on the canonical k-mer).
 * - Builds edges between minimizer-nodes if two k-mers (in the same read, after canonicalization)
 * would have been adjacent.
 * - Outputs unitigs (contiguous paths of minimizers, translated to sequences) in FASTA.
 *
 * Doesn’t implement tip‐removal, bubble‐popping, or abundance filtering.
 * This is a proof‐of‐concept skeleton, enhanced for correctness with canonical k-mers.
 */
public class MinimizerDBG {

    private final int k;           // k-mer length
    private final int w;           // window length for minimizer
    private final boolean verbose;

    // Maps each minimizer‐string (from a canonical k-mer) → its representative canonical k-mer (as String)
    private final Map<String, String> bucketedKmers = new HashMap<>();

    // Adjacency between minimizers: for each minimizer M, store set of neighbor minimizers N.
    private final Map<String, Set<String>> adj = new HashMap<>();

    public MinimizerDBG(int k, int w, boolean verbose) {
        if (w <= 0) {
            throw new IllegalArgumentException("Window size w must be positive.");
        }
        if (k <= 0) {
            throw new IllegalArgumentException("k-mer length k must be positive.");
        }
        if (w > k) {
            throw new IllegalArgumentException("Window size w cannot exceed k-mer length k.");
        }
        this.k = k;
        this.w = w;
        this.verbose = verbose;
    }

    // Helper function to get the reverse complement of a DNA sequence
    private String reverseComplement(String dna) {
        StringBuilder sb = new StringBuilder(dna.length());
        for (int i = dna.length() - 1; i >= 0; i--) {
            char base = dna.charAt(i);
            switch (base) {
                case 'A': sb.append('T'); break;
                case 'T': sb.append('A'); break;
                case 'C': sb.append('G'); break;
                case 'G': sb.append('C'); break;
                case 'N': sb.append('N'); break; // 'N' maps to 'N'
                default:
                    // Or throw an error for unexpected characters
                    sb.append(base); // Keep unknown characters as is for now
            }
        }
        return sb.toString();
    }

    // Helper function to get the canonical k-mer (lexicographically smaller of k-mer and its RC)
    private String getCanonical(String kmer) {
        String rcKmer = reverseComplement(kmer);
        return kmer.compareTo(rcKmer) <= 0 ? kmer : rcKmer;
    }

    /**
     * 1) Reads every FASTQ record.
     * 2) For each read, extracts all k-mers.
     * 3) For each k-mer, find its canonical form.
     * 4) Compute the minimizer of the canonical k-mer and add the canonical k-mer to the bucket.
     * 5) Record adjacency: if consecutive k-mers in the same read lead to different
     * (canonical) minimizers, link these minimizers.
     */
    public void buildFromFastq(File fastqFile) throws IOException {
        long readsProcessed = 0;
        long kmerCount = 0;

        try (FastqReader reader = new FastqReader(fastqFile)) {
            for (FastqRecord rec : reader) {
                readsProcessed++;
                String seq = rec.getReadString().toUpperCase();
                if (seq.length() < k) continue;

                String prevCanonMinimizer = null;

                for (int i = 0; i <= seq.length() - k; i++) {
                    String currentKmerOriginalStrand = seq.substring(i, i + k);

                    if (currentKmerOriginalStrand.contains("N")) {
                        prevCanonMinimizer = null; // Break chain of adjacency if 'N' is encountered
                        continue;
                    }
                    kmerCount++;

                    String canonicalCurrentKmer = getCanonical(currentKmerOriginalStrand);
                    String currentMinimizer = computeMinimizer(canonicalCurrentKmer);

                    // Bucket the canonical k-mer under its minimizer
                    bucketedKmers.putIfAbsent(currentMinimizer, canonicalCurrentKmer);

                    // Record adjacency between minimizers of canonical k-mers
                    if (prevCanonMinimizer != null && !prevCanonMinimizer.equals(currentMinimizer)) {
                        adj.computeIfAbsent(prevCanonMinimizer, key -> new HashSet<>()).add(currentMinimizer);
                        adj.computeIfAbsent(currentMinimizer, key -> new HashSet<>()).add(prevCanonMinimizer);
                    }
                    prevCanonMinimizer = currentMinimizer;
                }

                if (readsProcessed % 100000 == 0 && verbose) {
                    System.out.printf("  Processed %,d reads (%,d k-mers), minimizer buckets: %,d%n",
                            readsProcessed, kmerCount, bucketedKmers.size());
                }
            }
        }

        if (verbose) {
            System.out.printf("Finished processing %s.%nTotal reads: %,d, Total valid k-mers: %,d.%n",
                    fastqFile.getName(), readsProcessed, kmerCount);
            System.out.printf("Minimizer buckets: %,d, Adjacency entries: %,d.%n",
                    bucketedKmers.size(), adj.size());
        }
    }

    /** Find the lexicographically smallest substring of length w inside this (canonical) k-mer. */
    private String computeMinimizer(String canonicalKmer) {
        String bestMinimizer = null;
        // Iterate over all w-length windows in the k-mer
        for (int i = 0; i <= canonicalKmer.length() - w; i++) {
            String currentWindow = canonicalKmer.substring(i, i + w);
            if (bestMinimizer == null || currentWindow.compareTo(bestMinimizer) < 0) {
                bestMinimizer = currentWindow;
            }
        }
        return bestMinimizer;
    }

    /**
     * Write each “unitig” (maximal non-branching path of minimizers) to output FASTA.
     * Uses a greedy walk approach.
     */
    public void writeUnitigsAsFasta(File outFasta) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFasta))) {
            long contigCount = 0;
            Set<String> visitedMinimizers = new HashSet<>();

            // Process nodes that are potential path endpoints (degree != 2)
            for (String minimizer : bucketedKmers.keySet()) {
                if (visitedMinimizers.contains(minimizer)) continue;

                int degree = adj.getOrDefault(minimizer, Collections.emptySet()).size();
                if (degree == 2) continue; // Skip interior nodes for now; they'll be part of paths from endpoints

                List<String> path = walkFrom(minimizer, visitedMinimizers);
                if (!path.isEmpty()) {
                    String contigSeq = buildSequenceFromPath(path);
                    if (!contigSeq.isEmpty()) {
                        contigCount++;
                        writer.write(String.format(">unitig_%d_len_%d_path_%s%n%s%n",
                                contigCount, contigSeq.length(), pathToString(path), contigSeq));
                    }
                }
            }

            // Process remaining nodes, which must be part of cycles
            for (String minimizer : bucketedKmers.keySet()) {
                if (visitedMinimizers.contains(minimizer)) continue;
                // Any unvisited node here must be part of a cycle (all nodes in cycle have degree 2)
                List<String> cyclePath = walkCycle(minimizer, visitedMinimizers);
                if (!cyclePath.isEmpty()) {
                    String contigSeq = buildSequenceFromPath(cyclePath);
                    if (!contigSeq.isEmpty()) {
                        contigCount++;
                        writer.write(String.format(">unitig_%d_cycle_len_%d_path_%s%n%s%n",
                                contigCount, contigSeq.length(), pathToString(cyclePath), contigSeq));
                    }
                }
            }

            if (verbose) {
                System.out.printf("Wrote %,d unitigs to %s%n", contigCount, outFasta.getAbsolutePath());
            }
        }
    }

    private String pathToString(List<String> path) {
        if (path == null || path.isEmpty()) return "empty";
        return String.join("_", path); // Simple representation for debugging
    }


    /**
     * Starting from “startMinimizer,” walk in available directions through degree==2 nodes.
     * Marks all minimizers on the path as “visited.”
     */
    private List<String> walkFrom(String startMinimizer, Set<String> visited) {
        List<String> path = new ArrayList<>();
        path.add(startMinimizer);
        visited.add(startMinimizer);

        // Walk "forward" (one arbitrary direction from startMinimizer)
        String current = startMinimizer;
        String prev = null; // Stores the node we just came from to avoid immediate backtracking in chains
        while (true) {
            Set<String> neighbors = adj.getOrDefault(current, Collections.emptySet());
            if (adj.getOrDefault(current, Collections.emptySet()).size() != 2 && !current.equals(startMinimizer)) { // current is an endpoint
                break;
            }

            String next = null;
            int unvisitedNeighbors = 0;
            for (String neighbor : neighbors) {
                if (neighbor.equals(prev)) continue; // Don't go back immediately
                if (!visited.contains(neighbor)) {
                    next = neighbor;
                    unvisitedNeighbors++;
                } else if (neighbor.equals(startMinimizer) && path.size() > 1) { // Path loops back to start (not first step)
                    // This condition may not be strictly necessary if walkCycle handles all cycles
                    // but can terminate walks that close on themselves early.
                    next = null; // Force break
                    break;
                }
            }

            // If current is startMinimizer, it can have degree 1 (start of a simple path)
            // or degree > 2 (junction). We explore one path.
            // If current is an internal node (degree 2), it should have one unvisited neighbor (excluding prev).
            if (current.equals(startMinimizer)) {
                if (unvisitedNeighbors == 0) break; // No way forward from start
                // If multiple unvisited neighbors, pick one. (The loop above picks the last one found)
            } else { // current is an internal node of the path being built
                if (adj.getOrDefault(current, Collections.emptySet()).size() != 2) break; // Hit a branch or end
                if (unvisitedNeighbors != 1) break; // Should be exactly one way forward if degree 2 and not visited
            }


            if (next == null || visited.contains(next)) { // No valid next step or already visited
                break;
            }

            path.add(next);
            visited.add(next);
            prev = current;
            current = next;
        }


        // Walk “backward” from startMinimizer along other paths if startMinimizer is a junction
        // or to complete the other side if startMinimizer was an arbitrary internal point.
        // The current `walkFrom` is simplified to extend mainly in one "dominant" direction.
        // A more exhaustive bi-directional walk from general nodes might be needed for complex cases
        // or if startMinimizer is not a natural endpoint.
        // The original code had a more explicit two-directional walk from degree != 2 nodes.
        // For simplicity and matching the structure of finding paths from specific start types (endpoints/cycle entries):
        // This implementation starts a uni-directional walk. If the startMinimizer is degree 1, this is fine.
        // If it's a junction, it explores one branch. Other branches will be explored when their
        // turn comes in the main loop of writeUnitigsAsFasta.
        // The original backward walk code can be re-integrated if needed for a specific walk strategy.
        // The current path starts with startMinimizer and extends one way.

        // If startMinimizer had degree 1 and we walked from it, path is complete.
        // If startMinimizer was a junction (>2), this walk explores one branch.
        // The original code attempted a bi-directional walk from a single start point.
        // The following is closer to the original's second (backward) pass:
        current = startMinimizer;
        prev = path.size() > 1 ? path.get(1) : null; // Node already in path, "forward" from startMinimizer

        while(true) {
            Set<String> neighbors = adj.getOrDefault(current, Collections.emptySet());
            if (adj.getOrDefault(current, Collections.emptySet()).size() != 2 && !current.equals(startMinimizer)) {
                break;
            }

            String nextPrepend = null;
            int unvisitedPrependCandidates = 0;
            for(String neighbor : neighbors) {
                if (neighbor.equals(prev)) continue; // Don't go to the node we just "came from" in this direction
                if (!visited.contains(neighbor)) {
                    nextPrepend = neighbor;
                    unvisitedPrependCandidates++;
                }
            }

            if (current.equals(startMinimizer)) {
                // At startMinimizer, we look for another branch to prepend
                if (unvisitedPrependCandidates == 0) break;
            } else { // Internal node for the prepending part
                if (adj.getOrDefault(current, Collections.emptySet()).size() != 2) break;
                if (unvisitedPrependCandidates != 1) break; // Should be one way "backward"
            }


            if (nextPrepend == null || visited.contains(nextPrepend)) { // No valid node to prepend or already visited by another path
                break;
            }

            path.add(0, nextPrepend);
            visited.add(nextPrepend);
            prev = current;
            current = nextPrepend;
        }

        return path;
    }

    /**
     * Handle pure cycles: all nodes have degree==2.
     * Start from “minimizer,” walk until you return to it, marking visited.
     */
    private List<String> walkCycle(String startMinimizer, Set<String> visited) {
        List<String> cyclePath = new ArrayList<>();
        String current = startMinimizer;
        String prev = null; // To avoid immediate U-turns in simple 2-node cycles

        do {
            cyclePath.add(current);
            visited.add(current);

            Set<String> neighbors = adj.getOrDefault(current, Collections.emptySet());
            if (neighbors.size() != 2) { // Should not happen in a true cycle part
                System.err.println("Warning: Node " + current + " in supposed cycle has degree " + neighbors.size());
                // Potentially break or handle error, cycle might be malformed or path ended.
                return cyclePath; // Return what we have so far
            }

            String nextNode = null;
            for (String neighbor : neighbors) {
                if (!neighbor.equals(prev)) {
                    nextNode = neighbor;
                    break; // Pick the first valid neighbor to continue cycle
                }
            }

            if (nextNode == null) {
                // This could happen if prev was the only other node (e.g. 2-node cycle and picked wrong 'prev' initially)
                // Or if somehow both neighbors are 'prev' (logic error)
                // For a 2-node cycle (A-B-A), if current=A, prev=null, nextNode=B.
                // Then current=B, prev=A, nextNode should be A.
                if (neighbors.size() == 1 && prev != null) { // Only one neighbor and it's prev (dead end?)
                    System.err.println("Warning: Stuck in cycle walk at " + current);
                    return cyclePath;
                } else if (neighbors.size() > 0 && prev == null && cyclePath.size() > 1) { // Multiple neighbors but no valid 'next' found without prev
                    // This indicates an issue or that we've completed a loop back to a node that's not startMinimizer but connected
                }
                // If nextNode is still null, it means we can't move.
                if (nextNode == null && !neighbors.isEmpty()) { // Try picking any neighbor if stuck.
                    nextNode = neighbors.iterator().next();
                } else if (nextNode == null) {
                    System.err.println("Error: Cannot find next node in cycle from " + current);
                    return cyclePath; // Cannot continue
                }
            }
            prev = current;
            current = nextNode;

        } while (current != null && !current.equals(startMinimizer));

        // If cycle correctly closed, current should be startMinimizer. Add it to make the path explicit if needed by buildSequence.
        // However, buildSequenceFromPath often expects the sequence not to repeat the start/end kmer fully.
        // The current loop structure adds startMinimizer first, then proceeds until current becomes startMinimizer again, then exits.
        // So the path contains startMinimizer once at the beginning.

        return cyclePath;
    }


    /**
     * Given a list of minimizers [m₀, m₁, ..., mₚ], reconstruct a contig sequence.
     * Uses the representative canonical k-mer for m₀, then appends the last character
     * of the representative canonical k-mers for m₁, ..., mₚ.
     */
    private String buildSequenceFromPath(List<String> path) {
        if (path == null || path.isEmpty()) {
            return "";
        }

        String firstMinimizer = path.get(0);
        String representativeKmer = bucketedKmers.get(firstMinimizer);

        if (representativeKmer == null) {
            System.err.printf("Error: No k-mer found for the first minimizer in path: %s. Skipping path.%n", firstMinimizer);
            return "";
        }
        if (representativeKmer.length() != k) {
            System.err.printf("Error: Representative k-mer '%s' for minimizer '%s' has incorrect length %d (expected %d). Skipping path.%n",
                    representativeKmer, firstMinimizer, representativeKmer.length(), k);
            return "";
        }

        StringBuilder contig = new StringBuilder(representativeKmer);

        for (int i = 1; i < path.size(); i++) {
            String currentMinimizer = path.get(i);
            String nextRepresentativeKmer = bucketedKmers.get(currentMinimizer);

            if (nextRepresentativeKmer == null) {
                System.err.printf("Warning: No k-mer found for minimizer '%s' in path. Contig may be truncated.%n", currentMinimizer);
                break; // Stop extending this contig
            }
            if (nextRepresentativeKmer.length() != k) {
                System.err.printf("Warning: Representative k-mer '%s' for minimizer '%s' has incorrect length %d (expected %d). Contig may be truncated.%n",
                        nextRepresentativeKmer, currentMinimizer, nextRepresentativeKmer.length(), k);
                break; // Stop extending this contig
            }

            // Basic assumption: the representative k-mers are k-1 overlapping.
            // Append the last character of the next k-mer.
            contig.append(nextRepresentativeKmer.charAt(k - 1));
        }
        return contig.toString();
    }

    // Main method for testing (optional)
    public static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("Usage: java -cp <classpath> org.conncoll.MinimizerDBG <k> <w> <input.fastq> [output.fasta]");
            System.err.println("Example: java -cp . org.conncoll.MinimizerDBG 31 15 reads.fq unitigs.fa");
            System.exit(1);
        }

        try {
            int k = Integer.parseInt(args[0]);
            int w = Integer.parseInt(args[1]);
            File fastqFile = new File(args[2]);
            File outFile = new File(args.length > 3 ? args[3] : "unitigs_k" + k + "_w" + w + ".fasta");

            MinimizerDBG dbg = new MinimizerDBG(k, w, true); // verbose = true

            System.out.println("Building Minimizer de Bruijn Graph...");
            System.out.printf("Parameters: k=%d, w=%d, input FASTQ: %s%n", k, w, fastqFile.getAbsolutePath());

            long startTime = System.currentTimeMillis();
            dbg.buildFromFastq(fastqFile);
            long endTime = System.currentTimeMillis();
            System.out.printf("Graph construction took %.2f seconds.%n", (endTime - startTime) / 1000.0);

            System.out.println("Writing unitigs to FASTA: " + outFile.getAbsolutePath());
            startTime = System.currentTimeMillis();
            dbg.writeUnitigsAsFasta(outFile);
            endTime = System.currentTimeMillis();
            System.out.printf("Unitig writing took %.2f seconds.%n", (endTime - startTime) / 1000.0);

            System.out.println("Process completed.");

        } catch (NumberFormatException e) {
            System.err.println("Error: k and w must be integers.");
            e.printStackTrace();
        } catch (IllegalArgumentException e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
        }catch (IOException e) {
            System.err.println("File I/O error: " + e.getMessage());
            e.printStackTrace();
        } catch (Exception e) {
            System.err.println("An unexpected error occurred: " + e.getMessage());
            e.printStackTrace();
        }
    }
}