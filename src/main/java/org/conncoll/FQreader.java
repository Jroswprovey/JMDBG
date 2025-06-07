package org.conncoll;

import java.io.*;
import java.util.zip.GZIPInputStream;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqReader;

public class FQreader {

    private static File FQfile;

    public static void setFile(File InFile){
         FQfile = InFile;
    }

    public static void readRecord() throws IOException {

        try(FastqReader reader = new FastqReader(FQfile)){




        }




//        System.out.println("--- FASTQ: First Record Metrics ---");
//        System.out.println("Read Name: " + record.getReadString());
//        System.out.println("Read Length: " + record.length());
//        System.out.println("Read Bases (first 30): " + record.getReadString().substring(0, Math.min(30, record.length())));
//        System.out.println("Base Qualities (first 30): " + record.getBaseQualityString().substring(0, Math.min(30, record.getBaseQualityString().length())));

    }


}
