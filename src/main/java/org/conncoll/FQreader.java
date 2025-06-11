package org.conncoll;

import java.io.*;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqReader;

public class FQreader {

    private static File FQfile;

    public static void setFile(File InFile){
         FQfile = InFile;
    }

    public static String readRecord() throws IOException {

        try(FastqReader reader = new FastqReader(FQfile)){
           FastqRecord record = reader.next();
           return record.getReadString();
        }


    }


}
