package org.conncoll.MDBG;

/**
 * A utility class for handling 2-bit binary encoding of DNA k-mers.
 * This is a high-performance replacement for String-based k-mer operations.
 * A = 00, C = 01, G = 10, T = 11
 */
public final class KmerEncoder {

    // Make the constructor private so this utility class cannot be instantiated.
    private KmerEncoder() {}

    /**
     * Converts a DNA string k-mer into its 2-bit encoded long representation.
     * @param kmer The DNA string (e.g., "GAT").
     * @return A long representing the k-mer in binary.
     */
    public static long stringToLong(String kmer) {
        long binaryKmer = 0;
        for (char base : kmer.toCharArray()) {
            binaryKmer <<= 2; // Shift left by 2 to make room for the new base
            switch (base) {
                // 'A' is 00, so we do nothing (binaryKmer |= 0)
                case 'C' -> binaryKmer |= 1; // 01
                case 'G' -> binaryKmer |= 2; // 10
                case 'T' -> binaryKmer |= 3; // 11
            }
        }
        return binaryKmer;
    }

    /**
     * Converts a 2-bit encoded long back into a DNA string.
     * Essential for debugging and final output.
     * @param binaryKmer The long representation.
     * @param k The k-mer length.
     * @return The DNA string (e.g., "GAT").
     */
    public static String longToString(long binaryKmer, int k) {
        StringBuilder sb = new StringBuilder();
        long mask = 3; // A mask to extract the last 2 bits (binary 11)

        for (int i = 0; i < k; i++) {
            long twoBits = binaryKmer & mask; // Get the last 2 bits
            if (twoBits == 0) sb.append('A');
            else if (twoBits == 1) sb.append('C');
            else if (twoBits == 2) sb.append('G');
            else sb.append('T');
            binaryKmer >>= 2; // Shift right by 2 to process the next base
        }
        return sb.reverse().toString(); // The characters were added in reverse, so reverse the final string
    }

    /**
     * Calculates the reverse complement of a 2-bit encoded k-mer using bitwise operations.
     * @param binaryKmer The k-mer to complement.
     * @param k The length of the k-mer.
     * @return The reverse complement as a long.
     */
    public static long getReverseComplement(long binaryKmer, int k) {
        long revComp = 0;
        for (int i = 0; i < k; i++) {
            revComp <<= 2;
            long lastTwoBits = binaryKmer & 3;
            // Complementing (A<->T, C<->G) in 2-bit encoding is just flipping the bits (XOR with 11)
            revComp |= (lastTwoBits ^ 3);
            binaryKmer >>= 2;
        }
        return revComp;
    }

    /**
     * Finds the canonical representation of a 2-bit encoded k-mer.
     * This is beautifully simple: it's just the smaller of the two numbers.
     * @param binaryKmer The k-mer.
     * @param k The length of the k-mer.
     * @return The canonical k-mer as a long.
     */
    public static long getCanonical(long binaryKmer, int k) {
        long revComp = getReverseComplement(binaryKmer, k);
        return Math.min(binaryKmer, revComp);
    }
}