package org.conncoll.MDBG;

import java.util.*;


public class MinimizerDeBruijnGraph {

    private final int k;
    private final int w;
    private final Map<String, Set<String>> graph; // The final graph structure

    public MinimizerDeBruijnGraph(int k, int w){

        this.k = k;
        this.w = w;
        this.graph = new HashMap<>();

    }

    public List<String> getAllKmers(String sequence, int k){

        ArrayList<String> kmers = new ArrayList<>();

        for (int i = 0; i <= sequence.length() - k; i++) {

            kmers.add(sequence.substring(i,i + k));
        }
        return kmers;
    }



    private List<String> findMinimizers(List<String> kmers, int w){

        ArrayList<String> minimizers = new ArrayList<>();
        String curMin;

        for (int i = 0; i <= kmers.size() - w; i++) {
            curMin = kmers.get(i);

            for(int j = 0; j < kmers.subList(i,i + w).size(); j++){

                String canonicalWinner = DNAUtils.getLexicographicCanonical(curMin);
                String canonicalNext = DNAUtils.getLexicographicCanonical(kmers.get(i+j));

                if (canonicalNext.compareTo(canonicalWinner) < 0){

                    curMin = kmers.get(i+j);

                }
            }
            minimizers.add(curMin);
        }
        return minimizers;
    }



    private void buildGraph(List<String> minimizers){

    //Loop through all minimizes and links adjacent pairs
        for (int i = 0; i < minimizers.size()-1 ; i++) {
            String curNode = minimizers.get(i);
            String nextNode = minimizers.get(i+1);

            if (curNode.equals(nextNode)){
                continue;
            }


            this.graph.putIfAbsent(curNode, new HashSet<>());
            this.graph.get(curNode).add(nextNode);

        }
    }



    public void printGraph() {

        System.out.println("\n--- Minimizer de Bruijn Graph (Adjacency List) ---");
        if (this.graph.isEmpty()) {
            System.out.println("Graph is empty.");
            return;
        }
            // Loop through each node that has outgoing edges
        for (Map.Entry<String, Set<String>> entry : this.graph.entrySet()) {
            String sourceNode = entry.getKey();
            Set<String> neighbors = entry.getValue();

            // Print the node and its list of neighbors
            System.out.println(sourceNode + " -> " + neighbors);
        }
        System.out.println("----------------------------------------------------");
    }



    public String generateGraphvizDot() {

        StringBuilder dotBuilder = new StringBuilder();
        dotBuilder.append("digraph MDBG {\n");
        dotBuilder.append(" rankdir=LR;\n");
        dotBuilder.append(" node [shape=box, style=rounded];\n");

        // Loop through every edge in our graph

        for (Map.Entry<String, Set<String>> entry : this.graph.entrySet()) {

            String sourceNode = entry.getKey();

            for (String destNode : entry.getValue()) {

                // Add a line to the DOT string for each edge
                // e.g., " "GAT" -> "ATC";"

                dotBuilder.append(" \"" + sourceNode + "\" -> \"" + destNode + "\";\n");
            }
        }

        dotBuilder.append("}");
        return dotBuilder.toString();
    }



    public void addSequence(String sequence) {

        // These calls are for just this one sequence
        List<String> kmers = this.getAllKmers(sequence, this.k);
        List<String> minimizers = this.findMinimizers(kmers, this.w);
        this.buildGraph(minimizers); // buildGraph adds the new edges to the existing graph

    }



    private List<String> walkPath(String startNode, Set<String> visitedNodes){

        ArrayList<String> path = new ArrayList<>();

        String curNode = startNode;

        while(curNode != null){

            if(visitedNodes.contains(curNode)){
                break;
            }else {
                visitedNodes.add(curNode);
                path.add(curNode);
                Set<String> neighbors = this.graph.get(curNode);

                if (neighbors == null || neighbors.size() != 1){

                    curNode = null;

                }else {

                    curNode = neighbors.iterator().next();
                }
            }
        }
        return path;
    }





    public List<List<String>> assembleContigs() {

        List<List<String>> allContigs = new ArrayList<>();
        Set<String> visitedNodes = new HashSet<>();

        // We iterate through every node as a potential starting point

        for (String potentialStartNode : this.graph.keySet()) {
            if (visitedNodes.contains(potentialStartNode)) {
                continue; // Skip if we've already used this node
            }

            // If it's a new node, start walking a path
            List<String> newPath = walkPath(potentialStartNode, visitedNodes);
            if (!newPath.isEmpty()) {
                allContigs.add(newPath);
            }
        }
        return allContigs;
    }
}