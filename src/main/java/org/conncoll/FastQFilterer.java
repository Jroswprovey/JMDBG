package org.conncoll;

import htsjdk.samtools.fastq.*;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Set;
import java.util.zip.GZIPOutputStream;
import java.io.PrintStream;

public class FastQFilterer {


    public static File filterFastq(File inputFastqFile, File outputFastqFile, Set<String> readNamesToRemove) throws IOException {
        System.out.println("Starting" + inputFastqFile );

        try (FastqReader fastqReader = new FastqReader(inputFastqFile)) {

            OutputStream os;

            if(outputFastqFile.getName().endsWith(".gz")){
                os = new GZIPOutputStream(new FileOutputStream(outputFastqFile));
            }else {
                os = new FileOutputStream(outputFastqFile);
            }

            long totalReads;
            long writtenReads;
            FastqWriter basicWriter = new BasicFastqWriter(new PrintStream(os));
            try (FastqWriter fastqWriter = new AsyncFastqWriter(basicWriter, 500)) {


                totalReads = 0;
                writtenReads = 0;

                for (FastqRecord fastqRecord : fastqReader) {
                    totalReads++;

                    String currentFastqReadName = fastqRecord.getReadName();

                    if (!readNamesToRemove.contains(currentFastqReadName)) {
                        fastqWriter.write(fastqRecord);
                        writtenReads++;
                    }

                    if (totalReads % 1000000 == 0 && !Menu.verbose) {
                        System.out.println("Processed " + totalReads + " FASTQ records, written " + writtenReads + "...");
                    } else if (Menu.verbose) {
                        System.out.println("Processed " + totalReads + " FASTQ records, written " + writtenReads + "...");
                    }
                }
            }

            System.out.println("Finished processing FASTQ.");
            System.out.println("Total FASTQ records read: " + totalReads);
            System.out.println("Records written to filtered FASTQ: " + writtenReads);
            System.out.println("Records removed (found in SAM as mapped): " + (totalReads - writtenReads));

        }
        return outputFastqFile;
    }
}
