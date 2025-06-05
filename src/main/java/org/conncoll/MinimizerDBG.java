package org.conncoll;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.*;
import java.util.*;

/**
 * A very basic “minimizer de Bruijn graph” builder:
 *   - Buckets k-mers by their w-length minimizer.
 *   - Builds edges between minimizer-nodes if two k-mers (in the same read) overlap by k-1.
 *   - Outputs unitigs (contiguous paths) in FASTA.
 *
 * Doesn’t implement tip‐removal, bubble‐popping, or abundance filtering.
 * This is a proof‐of‐concept skeleton
 */
public class MinimizerDBG {

    private final int k;           // k-mer length
    private final int w;           // window length for minimizer
    private final boolean verbose;

    // Maps each minimizer‐string → a representative k-mer (as String) that belongs to that bucket
    private final Map<String, String> bucketedKmers = new HashMap<>();

    // Adjacency between minimizers: for each minimizer M, store set of neighbor minimizers N.
    private final Map<String, Set<String>> adj = new HashMap<>();

    public MinimizerDBG(int k, int w, boolean verbose) {
        if (w > k) {
            throw new IllegalArgumentException("Window size w cannot exceed k-mer length k");
        }
        this.k = k;
        this.w = w;
        this.verbose = verbose;
    }

    /**
     * 1) Reads every FASTQ record.
     * 2) For each read, extracts all k-mers.
     * 3) For each k-mer, find its minimizer and add to the bucket.
     * 4) Also record adjacency: consecutive k-mers in the same read → link their minimizers.
     */
    public void buildFromFastq(File fastqFile) throws IOException {
        long readsProcessed = 0;
        try (FastqReader reader = new FastqReader(fastqFile)) {
            for (FastqRecord rec : reader) {
                readsProcessed++;

                String seq = rec.getReadString().toUpperCase();
                if (seq.length() < k) continue; // skip too-short reads

                // Slide a window of length k across the read
                String prevMin = null;
                for (int i = 0; i <= seq.length() - k; i++) {
                    String kmer = seq.substring(i, i + k);
                    if (kmer.contains("N")) continue; // skip ambiguous bases

                    // 2a) find the minimizer of this k-mer
                    String minim = computeMinimizer(kmer);

                    // 2b) bucket the k-mer under that minimizer
                    bucketedKmers.putIfAbsent(minim, kmer);

                    // 2c) record adjacency (if there was a previous k-mer in this read)
                    if (prevMin != null && !prevMin.equals(minim)) {
                        // link prevMin ↔ minim
                        adj.computeIfAbsent(prevMin, x -> new HashSet<>()).add(minim);
                        adj.computeIfAbsent(minim, x -> new HashSet<>()).add(prevMin);
                    }

                    prevMin = minim;
                }

                if (readsProcessed % 100000 == 0 && verbose) {
                    System.out.printf("  Processed %,d reads, buckets: %,d%n", readsProcessed, bucketedKmers.size());
                }
            }
        }

        if (verbose) {
            System.out.printf("Finished bucketing k-mers: %d minimizer buckets, %d adjacency entries.%n",
                    bucketedKmers.size(), adj.size());
        }
    }

    /**
     * Write each “unitig” (maximal non-branching path of minimizers) to output FASTA.
     * Here, we do a simple “greedy walk” from each unvisited minimizer:
     *   1) If minimizer has degree ≠ 2, it’s a potential path endpoint.
     *   2) Walk from endpoints toward degree-2 nodes until you hit another endpoint.
     *   3) Convert that path of minimizers back to a nucleotide string by “splicing” their buckets:
     *      - We pick one representative k-mer from each bucket so that adjacent k-mers overlap by k-1.
     */
    public void writeUnitigsAsFasta(File outFasta) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFasta))) {
            long contigCount = 0;
            Set<String> visited = new HashSet<>();

            for (String minim : bucketedKmers.keySet()) {
                if (visited.contains(minim)) continue;

                // If degree != 2, treat this minimizer as a path “start” (or a single-node path).
                int deg = adj.getOrDefault(minim, Collections.emptySet()).size();
                if (deg == 2) continue; // skip interior nodes—will be handled from some endpoint

                // Now walk forward from “minim”:
                List<String> path = walkFrom(minim, visited);

                // Convert that path of minimizers → contig sequence
                String contigSeq = buildSequenceFromPath(path);
                contigCount++;
                writer.write(String.format(">contig_%d_len_%d%n%s%n",
                        contigCount, contigSeq.length(), contigSeq));
            }

            // Check for any cycles (all nodes in a cycle have degree=2 and won't be reached above).
            for (String minim : bucketedKmers.keySet()) {
                if (visited.contains(minim)) continue;
                // This is a cycle. Just pick any node in the cycle to walk around until we return.
                List<String> cyclePath = walkCycle(minim, visited);
                String contigSeq = buildSequenceFromPath(cyclePath);
                contigCount++;
                writer.write(String.format(">contig_%d_cycle_len_%d%n%s%n",
                        contigCount, contigSeq.length(), contigSeq));
            }

            if (verbose) {
                System.out.printf("Wrote %,d unitigs to %s%n", contigCount, outFasta);
            }
        }
    }

    /** Find the lexicographically smallest substring of length w inside this k-mer. */
    private String computeMinimizer(String kmer) {
        String best = null;
        for (int i = 0; i <= k - w; i++) {
            String sub = kmer.substring(i, i + w);
            if (best == null || sub.compareTo(best) < 0) {
                best = sub;
            }
        }
        return best;
    }

    /**
     * Starting from “startMinim,” walk in both directions through degree==2 nodes until endpoints.
     * Mark everything on that path “visited.”
     */
    private List<String> walkFrom(String startMinim, Set<String> visited) {
        List<String> path = new ArrayList<>();
        path.add(startMinim);
        visited.add(startMinim);

        // Walk “forward” (arbitrary orientation) until you hit a node with deg≠2 or a visited node.
        String curr = startMinim;
        String prev = null;
        while (true) {
            Set<String> neighbors = adj.getOrDefault(curr, Collections.emptySet());
            if (neighbors.size() != 2) break;

            // choose the neighbor that isn’t “prev”
            String next = null;
            for (String n : neighbors) {
                if (!n.equals(prev)) {
                    next = n;
                    break;
                }
            }
            if (next == null || visited.contains(next)) break;
            path.add(next);
            visited.add(next);
            prev = curr;
            curr = next;
        }

        // Walk “backward” from startMinim along the other side
        curr = startMinim;
        prev = null;
        // But now we look at neighbors and pick the one not already in “path” (if deg==2)
        while (true) {
            Set<String> neighbors = adj.getOrDefault(curr, Collections.emptySet());
            if (neighbors.size() != 2) break;

            String other = null;
            for (String n : neighbors) {
                if (!n.equals(prev) && !visited.contains(n)) {
                    other = n;
                    break;
                }
            }
            if (other == null) break;
            // prepend to path
            path.add(0, other);
            visited.add(other);
            prev = curr;
            curr = other;
        }

        return path;
    }

    /**
     * Handle pure cycles: all nodes have degree==2, so they never appear as a “start” above.
     * Start from “minim,” walk until you return to it.
     */
    private List<String> walkCycle(String minim, Set<String> visited) {
        List<String> cycle = new ArrayList<>();
        String start = minim;
        String curr = minim;
        String prev = null;
        do {
            cycle.add(curr);
            visited.add(curr);

            // pick the next neighbor (any one not equal to prev)
            String next = null;
            for (String n : adj.getOrDefault(curr, Collections.emptySet())) {
                if (!n.equals(prev)) {
                    next = n;
                    break;
                }
            }
            prev = curr;
            curr = next;
        } while (curr != null && !curr.equals(start));

        return cycle;
    }

    /**
     * Given a list of minimizers [m₀, m₁, ..., mₙ], reconstruct a contig sequence by
     * picking one k-mer from bucket[m₀], then extending with overlapping k-mers from each bucket in order.
     * Here we simply pick _any_ k-mer from each bucket that has the correct prefix/suffix match.
     */
    private String buildSequenceFromPath(List<String> path) {
        if (path.isEmpty()) return "";

        // 1) pick a representative k-mer from the first bucket:
        String firstMinim = path.get(0);
        String contig = bucketedKmers.get(firstMinim);

        // 2) now for each subsequent minimizer mᵢ, pick a k-mer from bucket that overlaps contig’s last k-1
        for (int i = 1; i < path.size(); i++) {
            String currMin = path.get(i);
            String candidate = bucketedKmers.get(currMin);
            String suffix = contig.substring(contig.length() - (k - 1)); // last k-1 bases
            if (candidate.startsWith(suffix)) {
                contig += candidate.charAt(k - 1);
            } else {
                contig += candidate.substring(k - 1);
            }
        }

        return contig;
    }
}