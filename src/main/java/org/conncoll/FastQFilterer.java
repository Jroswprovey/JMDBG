package org.conncoll;

import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

import java.io.File;
import java.io.IOException;
import java.util.Set;

public class FastQFilterer {


    public static void filterFastq(File inputFastqFile, File outputFastqFile, Set<String> readNamesToRemove) throws IOException {
        System.out.println("Starting" + inputFastqFile );
        FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
        try (FastqReader fastqReader = new FastqReader(inputFastqFile);
             FastqWriter fastqWriter = fastqWriterFactory.newWriter(outputFastqFile)) {

            long totalReads = 0;
            long writtenReads = 0;

            for (FastqRecord fastqRecord : fastqReader) {
                totalReads++;

                String currentFastqReadName = fastqRecord.getReadName();

                if (!readNamesToRemove.contains(currentFastqReadName)) {
                    fastqWriter.write(fastqRecord);
                    writtenReads++;
                }

                if (totalReads % 1000000 == 0 && !Menu.verbose) {
                    System.out.println("Processed " + totalReads + " FASTQ records, written " + writtenReads + "...");
                } else if (Menu.verbose){
                    System.out.println("Processed " + totalReads + " FASTQ records, written " + writtenReads + "...");
                }
            }

            System.out.println("Finished processing FASTQ.");
            System.out.println("Total FASTQ records read: " + totalReads);
            System.out.println("Records written to filtered FASTQ: " + writtenReads);
            System.out.println("Records removed (found in SAM as mapped): " + (totalReads - writtenReads));
        }
    }
}
