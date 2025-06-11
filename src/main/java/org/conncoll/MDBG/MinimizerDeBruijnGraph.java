package org.conncoll.MDBG;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.*;

public class MinimizerDeBruijnGraph {

    private final int k;
    private final int w;
    private final Map<Long, List<Edge>> graph;

    // The Edge record now stores the destination node as a long
    public record Edge(long destinationNode, long offset, int length) {}

    public MinimizerDeBruijnGraph(int k, int w) {
        this.k = k;
        this.w = w;
        this.graph = new HashMap<>();
    }

    // This is now an internal record to hold a minimizer and its original position
    private record MinimizerInfo(long sequence, int position) {}

    public void addSequence(String sequence, RandomAccessFile edgeFile) throws IOException {
        if (sequence.length() < k) return;

        List<MinimizerInfo> minimizers = this.findMinimizers(sequence, k, w);
        this.buildGraph(minimizers, sequence, edgeFile);
    }

    private List<MinimizerInfo> findMinimizers(String sequence, int k, int w) {
        ArrayList<MinimizerInfo> minimizers = new ArrayList<>();
        List<Long> kmers = new ArrayList<>();
        List<Long> canonicalKmers = new ArrayList<>();

        // Prime the k-mer list
        for (int i = 0; i <= sequence.length() - k; i++) {
            long kmerAsLong = KmerEncoder.stringToLong(sequence.substring(i, i + k));
            kmers.add(kmerAsLong);
            canonicalKmers.add(KmerEncoder.getCanonical(kmerAsLong, k));
        }

        // Find minimizers
        for (int i = 0; i <= kmers.size() - w; i++) {
            long windowMinimizer = -1;
            int windowMinimizerPosition = -1;

            for (int j = 0; j < w; j++) {
                long currentCanonical = canonicalKmers.get(i + j);
                if (windowMinimizer == -1 || currentCanonical < windowMinimizer) {
                    windowMinimizer = currentCanonical;
                    windowMinimizerPosition = i + j;
                }
            }
            minimizers.add(new MinimizerInfo(kmers.get(windowMinimizerPosition), windowMinimizerPosition));
        }
        return minimizers;
    }

    private void buildGraph(List<MinimizerInfo> minimizers, String originalRead, RandomAccessFile edgeFile) throws IOException {
        if (minimizers.size() < 2) {
            return;
        }

        for (int i = 0; i < minimizers.size() - 1; i++) {
            MinimizerInfo sourceInfo = minimizers.get(i);
            MinimizerInfo destInfo = minimizers.get(i + 1);

            long sourceNode = KmerEncoder.getCanonical(sourceInfo.sequence(), k);
            long destNode = KmerEncoder.getCanonical(destInfo.sequence(), k);

            if (sourceNode == destNode) {
                continue;
            }

            int pos1 = sourceInfo.position();
            int pos2 = destInfo.position();

            if (pos2 > pos1) { // Ensure they are in the correct order from the read
                String edgeSequence = originalRead.substring(pos1, pos2 + this.k);
                long offset = edgeFile.length();
                edgeFile.write(edgeSequence.getBytes());
                int length = edgeSequence.length();

                Edge pointerEdge = new Edge(destNode, offset, length);
                this.graph.putIfAbsent(sourceNode, new ArrayList<>());
                this.graph.get(sourceNode).add(pointerEdge);
            }
        }
    }

    private List<Long> walkPath(long startNode, Set<Long> visitedNodes) {
        ArrayList<Long> path = new ArrayList<>();
        Long curNode = startNode;

        while (curNode != null) {
            if (visitedNodes.contains(curNode)) {
                break;
            }
            visitedNodes.add(curNode);
            path.add(curNode);

            List<Edge> edges = this.graph.get(curNode);
            if (edges == null || edges.size() != 1) {
                curNode = null;
            } else {
                curNode = edges.get(0).destinationNode();
            }
        }
        return path;
    }

    public List<List<Long>> assembleContigs() {
        List<List<Long>> allContigs = new ArrayList<>();
        Set<Long> visitedNodes = new HashSet<>();

        for (long potentialStartNode : this.graph.keySet()) {
            if (visitedNodes.contains(potentialStartNode)) {
                continue;
            }
            List<Long> newPath = walkPath(potentialStartNode, visitedNodes);
            if (!newPath.isEmpty()) {
                allContigs.add(newPath);
            }
        }
        return allContigs;
    }

    public String stitchContig(List<Long> contigNodePath, File edgeFile) throws IOException {
        if (contigNodePath == null || contigNodePath.size() < 2) {
            return contigNodePath != null && !contigNodePath.isEmpty() ? KmerEncoder.longToString(contigNodePath.get(0), k) : "";
        }

        StringBuilder finalSequence = new StringBuilder();

        try (RandomAccessFile reader = new RandomAccessFile(edgeFile, "r")) {
            for (int i = 0; i < contigNodePath.size() - 1; i++) {
                long sourceNode = contigNodePath.get(i);
                long destNode = contigNodePath.get(i + 1);

                List<Edge> edges = this.graph.get(sourceNode);
                if (edges == null) continue;

                for (Edge edge : edges) {
                    if (edge.destinationNode() == destNode) {
                        reader.seek(edge.offset());
                        byte[] buffer = new byte[edge.length()];
                        reader.readFully(buffer);
                        String edgeSequence = new String(buffer);

                        if (finalSequence.length() == 0) {
                            finalSequence.append(edgeSequence);
                        } else {
                            finalSequence.append(edgeSequence.substring(this.k));
                        }
                        break;
                    }
                }
            }
        }
        return finalSequence.toString();
    }

    public void printGraph() {
        System.out.println("\n--- Minimizer de Bruijn Graph (Adjacency List) ---");
        if (this.graph.isEmpty()) {
            System.out.println("Graph is empty.");
            return;
        }
        for (Map.Entry<Long, List<Edge>> entry : this.graph.entrySet()) {
            String sourceString = KmerEncoder.longToString(entry.getKey(), k);
            System.out.println(sourceString + " -> " + entry.getValue());
        }
        System.out.println("----------------------------------------------------");
    }

    public String generateGraphvizDot() {
        StringBuilder dotBuilder = new StringBuilder();
        dotBuilder.append("digraph MDBG {\n");
        dotBuilder.append("    rankdir=LR;\n");
        dotBuilder.append("    node [shape=box, style=rounded];\n");

        for (Map.Entry<Long, List<Edge>> entry : this.graph.entrySet()) {
            String sourceNode = KmerEncoder.longToString(entry.getKey(), k);
            for (Edge edge : entry.getValue()) {
                String destNode = KmerEncoder.longToString(edge.destinationNode(), k);
                String label = edge.length() + "bp";
                dotBuilder.append("    \"" + sourceNode + "\" -> \"" + destNode + "\" [label=\"" + label + "\"];\n");
            }
        }
        dotBuilder.append("}");
        return dotBuilder.toString();
    }
}