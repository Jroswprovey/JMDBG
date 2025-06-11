package org.conncoll.MDBG;

public class DNAUtils {

    public static char[] COMPLIMENTS = new char[128];

    static {
        //Dna compliment bases lookup
        COMPLIMENTS['A'] = 'T';
        COMPLIMENTS['T'] = 'A';
        COMPLIMENTS['G'] = 'C';
        COMPLIMENTS['C'] = 'G';
    }


    private static String getReverseComplement(String kmer){
        StringBuilder sb = new StringBuilder(kmer);
        StringBuilder rc = new StringBuilder();
        sb.reverse();

        for(int i = 0; i<sb.length(); ++i){
            rc.append(COMPLIMENTS[sb.charAt(i)]);
        }

        return rc.toString();
    }

    public static String getLexicographicCanonical(String kmer){
        String rcKmer = getReverseComplement(kmer);

        if (rcKmer.compareTo(kmer) < 0) {
            return rcKmer;
        } else {
            return kmer;
        }

    }

}
