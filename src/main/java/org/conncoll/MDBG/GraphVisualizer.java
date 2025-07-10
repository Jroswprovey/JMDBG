package org.conncoll.MDBG;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// Add this class to your project or put the method inside your existing class.
public class GraphVisualizer {

    /**
     * Prints a human-readable representation of the de Bruijn graph.
     *
     * @param graph         The graph structure (Map<Integer, List<Integer>>).
     * @param minimizerToId The map from k-mer strings to their integer IDs.
     */
    public static void printGraph(Map<Integer, List<Integer>> graph, Map<String, Integer> minimizerToId) {
        System.out.println("\n--- De Bruijn Graph Visualization ---");

        // Step 1: Create a reverse map for easy lookup of k-mer strings from IDs.
        // This is much more efficient than searching the original map every time.
        Map<Integer, String> idToMinimizer = new HashMap<>();
        for (Map.Entry<String, Integer> entry : minimizerToId.entrySet()) {
            idToMinimizer.put(entry.getValue(), entry.getKey());
        }

        // Check if the graph is empty
        if (graph.isEmpty()) {
            System.out.println("Graph is empty. No connections were made.");
            System.out.println("-------------------------------------\n");
            return;
        }

        // To make the output tidy and predictable, let's sort the nodes by their ID
        List<Integer> sortedNodeIds = new ArrayList<>(graph.keySet());
        Collections.sort(sortedNodeIds);

        // Step 2: Iterate through each node in the graph that has outgoing edges.
        for (Integer fromNodeId : sortedNodeIds) {
            String fromKmer = idToMinimizer.get(fromNodeId);

            // Build the output string for the current node
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("Node %d (%s) --> ", fromNodeId, fromKmer));

            // Step 3: Get its list of neighbors (the "to" nodes).
            List<Integer> neighbors = graph.get(fromNodeId);

            List<String> neighborInfo = new ArrayList<>();
            for (Integer toNodeId : neighbors) {
                String toKmer = idToMinimizer.get(toNodeId);
                neighborInfo.add(String.format("%d (%s)", toNodeId, toKmer));
            }
            // Join all the neighbor info with commas for a clean look
            sb.append(String.join(", ", neighborInfo));

            // Step 4: Print the final formatted line for this node.
            System.out.println(sb.toString());
        }
        System.out.println("-------------------------------------\n");
    }
}
