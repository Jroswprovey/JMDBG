package org.conncoll.MDBG;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Assembly {

    // A simple record to hold edge information.
    // This could also be a private static class inside Assembly.
    record Edge(int toId, String sequence) {}

    /**
     * Assembles unitigs by streaming from a sorted edge file, using pre-calculated node degrees.
     * This approach is memory-efficient as it does not load the entire graph into RAM.
     *
     * @param sortedEdgesPath The path to the file containing sorted graph edges.
     * @param outputPath      The path for the output FASTA file of contigs.
     * @param kmerSize        The k-mer size used during graph construction.
     * @param inDegrees       A map of node ID to its in-degree.
     * @param outDegrees      A map of node ID to its out-degree.
     * @throws IOException If an I/O error occurs.
     */
    public void assembleStreaming(String sortedEdgesPath, String outputPath, int kmerSize,
                                  Int2IntOpenHashMap inDegrees, Int2IntOpenHashMap outDegrees) throws IOException {

        System.out.println("Assembling unitigs via streaming...");

        // Key: The ID of the LAST node in a path. Value: The sequence of the path so far.
        Map<Integer, StringBuilder> activePaths = new HashMap<>();
        List<String> finalUnitigs = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(sortedEdgesPath))) {
            String line;
            int lastFromId = -1;
            List<Edge> currentBlock = new ArrayList<>(); // Holds all edges for the current fromId

            while ((line = reader.readLine()) != null) {
                String[] parts = line.split("\t");
                int fromId = Integer.parseInt(parts[0]);

                if (fromId != lastFromId && lastFromId != -1) {
                    // We've finished reading a block of edges for a node, so process it.
                    processBlock(lastFromId, currentBlock, activePaths, finalUnitigs, inDegrees, outDegrees, kmerSize);
                    currentBlock.clear();
                }

                currentBlock.add(new Edge(Integer.parseInt(parts[1]), parts[2]));
                lastFromId = fromId;
            }
            // Process the very last block in the file after the loop finishes.
            if (!currentBlock.isEmpty()) {
                processBlock(lastFromId, currentBlock, activePaths, finalUnitigs, inDegrees, outDegrees, kmerSize);
            }
        }

        // Any paths remaining in activePaths are also complete unitigs.
        for (StringBuilder sb : activePaths.values()) {
            finalUnitigs.add(sb.toString());
        }

        // Write all finalized unitigs to the output FASTA file.
        writeFasta(outputPath, finalUnitigs);
        System.out.println("Assembly complete. Wrote " + finalUnitigs.size() + " unitigs.");
    }

    /**
     * Processes a block of edges all starting from the same node.
     */
    private void processBlock(int fromId, List<Edge> edges, Map<Integer, StringBuilder> activePaths,
                              List<String> finalUnitigs, Int2IntOpenHashMap inDegrees, Int2IntOpenHashMap outDegrees, int kmerSize) {

        // A path that was being extended has now arrived at this node block.
        // Because this node is complex (a branch or merge point), the incoming path is now complete.
        StringBuilder incomingPath = activePaths.remove(fromId);
        if (incomingPath != null) {
            finalUnitigs.add(incomingPath.toString());
        }

        // Now, check if this node is the START of new unitigs.
        // This happens if it's a "head" node (not a simple 1-in, 1-out connector).
        // The simplest definition of a head is a node that doesn't have a simple "in" path.
        if (inDegrees.get(fromId) != 1) {
            // Start a new unitig for each outgoing edge.
            for (Edge edge : edges) {
                extendPath(new StringBuilder(edge.sequence), edge, activePaths, finalUnitigs, inDegrees, outDegrees, kmerSize);
            }
        }
    }

    /**
     * Extends a given path with a new edge and decides if the path is complete or should continue.
     */
    private void extendPath(StringBuilder currentSeq, Edge edge, Map<Integer, StringBuilder> activePaths,
                            List<String> finalUnitigs, Int2IntOpenHashMap inDegrees, Int2IntOpenHashMap outDegrees, int kmerSize) {

        // Append the non-overlapping part of the new edge's sequence.
        currentSeq.append(edge.sequence.substring(kmerSize));
        int toId = edge.toId;

        // Decide if the path is finished or should be put back in the active map.
        if (outDegrees.get(toId) == 1 && inDegrees.get(toId) == 1) {
            // This is a simple path, so put it back in activePaths to be extended later.
            activePaths.put(toId, currentSeq);
        } else {
            // This path hits a branch or a dead end, so it's a final unitig.
            finalUnitigs.add(currentSeq.toString());
        }
    }

    /**
     * Writes a list of sequences to a FASTA formatted file.
     */
    private void writeFasta(String outputPath, List<String> unitigs) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {
            for (int i = 0; i < unitigs.size(); i++) {
                String sequence = unitigs.get(i);
                writer.write(">unitig_" + i + " length_" + sequence.length() + "\n");

                // Wrap lines for FASTA format (e.g., every 80 characters)
                for (int j = 0; j < sequence.length(); j += 80) {
                    writer.write(sequence.substring(j, Math.min(j + 80, sequence.length())));
                    writer.newLine();
                }
            }
        }
    }
}