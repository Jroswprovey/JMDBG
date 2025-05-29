package org.conncoll;

import java.io.*;
import java.util.Spliterator;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import java.util.zip.GZIPInputStream;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqReader;

public class FQreader {

    private static File FQfile;

    public static void setFile(File InFile){
         FQfile = InFile;
    }

    public static void read() throws IOException {

        FileInputStream fis = null;
        try {
            fis = new FileInputStream(FQfile);
        } catch (FileNotFoundException e) {
            System.err.println("File not found");
        }

        if (fis == null) throw new AssertionError();

        GZIPInputStream gis = new GZIPInputStream(fis);
        InputStreamReader isr = new InputStreamReader(gis);
        BufferedReader br = new BufferedReader(isr);
        FastqReader fr = new FastqReader(FQfile,br);

        FastqRecord record = fr.next();



        System.out.println("--- FASTQ: First Record Metrics ---");
        System.out.println("Read Name: " + record.getReadString());
        System.out.println("Read Length: " + record.length());
        System.out.println("Read Bases (first 30): " + record.getReadString().substring(0, Math.min(30, record.length())));
        System.out.println("Base Qualities (first 30): " + record.getBaseQualityString().substring(0, Math.min(30, record.getBaseQualityString().length())));
    }


}
