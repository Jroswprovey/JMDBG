package org.conncoll;

import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class FastQFilterer {

    public static Set<String> getMappedReadNames(File samFile) throws IOException {
        Set<String> mappedReadNames = new HashSet<>();
        try (SamReader samReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(samFile)){

            for(final SAMRecord samRecord : samReader){
                //Loops through the records and if they are mapped, add them to the hashset
                if (!samRecord.getReadUnmappedFlag()){
                    mappedReadNames.add(samRecord.getReadName());
                }
            }
        }

        System.out.println("Found " + mappedReadNames.size() + " Mapped Reads");
        return mappedReadNames;

    }

    public static void filterFastq(File inputFastqFile, File outputFastqFile, Set<String> readNamesToRemove) throws IOException {
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

                if (totalReads % 1000000 == 0) {
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
