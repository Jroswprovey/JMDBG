package org.conncoll.MDBG;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import it.unimi.dsi.fastutil.longs.LongSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Assembly {

    record Edge(int toId, String sequence) {}
    record BubblePath(String sequence, int endNodeId) {}

    private static final int MAX_BUBBLE_DEPTH = 20;
    private static final int MAX_BUBBLE_LENGTH = 1500;

    public void assembleStreaming(String sortedEdgesPath, String outputPath, int kmerSize,
                                  Int2IntOpenHashMap inDegrees, Int2IntOpenHashMap outDegrees) throws IOException {

        System.out.println("Building in-memory graph...");
        Map<Integer, List<Edge>> graph = loadGraph(sortedEdgesPath);
        System.out.println("Graph loaded with " + graph.size() + " nodes that have outgoing edges.");

        LongSet traversedEdges = new LongOpenHashSet();
        List<String> finalUnitigs = new ArrayList<>();

        System.out.println("Assembling unitigs...");
        IntSet allNodes = new IntOpenHashSet(inDegrees.keySet());
        allNodes.addAll(outDegrees.keySet());

        for (int nodeId : allNodes) {
            // A unitig starts at a node that is NOT a simple corridor.
            if (inDegrees.get(nodeId) != 1 || outDegrees.get(nodeId) != 1) {
                List<Edge> outgoingEdges = graph.get(nodeId);
                if (outgoingEdges != null) {
                    for (Edge edge : outgoingEdges) {
                        walkPath(nodeId, edge, graph, traversedEdges, finalUnitigs, inDegrees, outDegrees, kmerSize);
                    }
                }
            }
        }

        writeFasta(outputPath, finalUnitigs, kmerSize);
        System.out.println("Assembly complete. Wrote " + finalUnitigs.size() + " unitigs.");
    }

    private long getEdgeId(int from, int to) {
        return ((long) from << 32) | (to & 0xFFFFFFFFL);
    }

    private void walkPath(int startNode, Edge firstEdge, Map<Integer, List<Edge>> graph, LongSet traversedEdges,
                          List<String> finalUnitigs, Int2IntOpenHashMap inDegrees, Int2IntOpenHashMap outDegrees, int kmerSize) {

        long firstEdgeId = getEdgeId(startNode, firstEdge.toId);
        if (traversedEdges.contains(firstEdgeId)) {
            return;
        }

        StringBuilder currentSequence = new StringBuilder(firstEdge.sequence);
        traversedEdges.add(firstEdgeId);

        int currentNodeId = firstEdge.toId;

        while (inDegrees.get(currentNodeId) == 1 && outDegrees.get(currentNodeId) == 1) {
            List<Edge> nextEdges = graph.get(currentNodeId);
            if (nextEdges == null || nextEdges.isEmpty()) break;

            Edge nextEdge = nextEdges.get(0);
            long nextEdgeId = getEdgeId(currentNodeId, nextEdge.toId);

            if (traversedEdges.contains(nextEdgeId)) {
                break;
            }

            // The sequence on the edge is the string connecting two minimizers.
            // We need to append the part of this string that comes after the overlapping k-mer.
            currentSequence.append(nextEdge.sequence.substring(kmerSize));

            traversedEdges.add(nextEdgeId);
            currentNodeId = nextEdge.toId;
        }

        // Bubble popping logic can be added back here later if needed.

        if (currentSequence.length() >= kmerSize) {
            finalUnitigs.add(currentSequence.toString());
        }
    }

    private BubblePath resolveBubble(int startNodeId, Map<Integer, List<Edge>> graph, LongSet traversedEdges,
                                     Int2IntOpenHashMap inDegrees, Int2IntOpenHashMap outDegrees, int kmerSize) {

        List<Edge> branches = graph.get(startNodeId);
        if (branches == null || branches.size() < 2) return null;

        List<BubblePath> resolvedPaths = new ArrayList<>();
        int consensusEndNode = -1;

        for (Edge branchEdge : branches) {
            long branchEdgeId = getEdgeId(startNodeId, branchEdge.toId);
            if (traversedEdges.contains(branchEdgeId)) continue;

            StringBuilder pathSequence = new StringBuilder(branchEdge.sequence);
            int pathCurrentNodeId = branchEdge.toId;
            boolean pathIsValid = true;

            for (int depth = 0; depth < MAX_BUBBLE_DEPTH; depth++) {
                if (pathSequence.length() > MAX_BUBBLE_LENGTH) {
                    pathIsValid = false;
                    break;
                }
                if (inDegrees.get(pathCurrentNodeId) > 1) break;
                if (outDegrees.get(pathCurrentNodeId) != 1) {
                    pathIsValid = false;
                    break;
                }

                List<Edge> nextEdges = graph.get(pathCurrentNodeId);
                if (nextEdges == null || nextEdges.isEmpty()) {
                    pathIsValid = false;
                    break;
                }
                Edge nextEdge = nextEdges.get(0);

                if (traversedEdges.contains(getEdgeId(pathCurrentNodeId, nextEdge.toId()))) {
                    pathIsValid = false;
                    break;
                }
                pathSequence.append(nextEdge.sequence.substring(kmerSize));
                pathCurrentNodeId = nextEdge.toId;
            }
            if (pathIsValid) {
                if (consensusEndNode == -1) consensusEndNode = pathCurrentNodeId;
                else if (consensusEndNode != pathCurrentNodeId) return null;

                resolvedPaths.add(new BubblePath(pathSequence.toString(), pathCurrentNodeId));
            } else {
                return null;
            }
        }
        if (resolvedPaths.size() != branches.size() || resolvedPaths.isEmpty()) return null;

        resolvedPaths.sort(Comparator.comparing(p -> p.sequence));
        BubblePath bestPath = resolvedPaths.get(0);

        return new BubblePath(bestPath.sequence.substring(kmerSize), bestPath.endNodeId);
    }

    private Map<Integer, List<Edge>> loadGraph(String sortedEdgesPath) throws IOException {
        Map<Integer, List<Edge>> graph = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(sortedEdgesPath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split("\t");
                int fromId = Integer.parseInt(parts[0]);
                int toId = Integer.parseInt(parts[1]);
                String sequence = parts[2];
                graph.computeIfAbsent(fromId, k -> new ArrayList<>()).add(new Edge(toId, sequence));
            }
        }
        return graph;
    }

    private void writeFasta(String outputPath, List<String> unitigs, int kmerSize) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {
            int writtenCount = 0;
            for (String sequence : unitigs) {
                if (sequence.length() < kmerSize) continue;
                writer.write(">unitig_" + writtenCount + " length_" + sequence.length() + "\n");
                for (int j = 0; j < sequence.length(); j += 80) {
                    writer.write(sequence.substring(j, Math.min(j + 80, sequence.length())));
                    writer.newLine();
                }
                writtenCount++;
            }
        }
    }
}