package org.conncoll;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;


public class samUtilities {

    public static Set<String> getMappedReadNames(File samFile) throws IOException {
        System.out.println("Getting mapped values from sam File");
        Set<String> mappedReadNames = new HashSet<>();
        try (SamReader samReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(samFile)) {

            for (final SAMRecord samRecord : samReader) {
                //Loops through the records and if they are mapped, add them to the hashset
                if (!samRecord.getReadUnmappedFlag()) {
                    mappedReadNames.add(samRecord.getReadName());
                }
            }
        }

        System.out.println("Found " + mappedReadNames.size() + " Mapped Reads");
        return mappedReadNames;

    }


    /*
     * getMaxMapQ will find the maximum MAPQ quality so
     * the user can be given some context when selecting a threshold
     */
    public static int getMaxMapQ(SAMRecordIterator recordIterator){
        int max = 0;

        while(recordIterator.hasNext()){
            SAMRecord record = recordIterator.next();
            if(record.getMappingQuality() > max){
             max = record.getMappingQuality();
            }
        }
        return max;
    }


    public static boolean isFileSorted(File samOrBamFile) throws IOException {

        if(!samOrBamFile.exists()){
            System.err.println("Error: File not found at: " + samOrBamFile.getAbsolutePath());
            throw new RuntimeException();
        }

        try(SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(samOrBamFile)) {

            //Defines the header then retrieves the sort order defined in the header
            SAMFileHeader header = reader.getFileHeader();
            SAMFileHeader.SortOrder sortOrder = header.getSortOrder();

            if(sortOrder == SAMFileHeader.SortOrder.unsorted){
                System.out.println("File is unsorted");
                return true; //This returns true because if its unsorted then it's technically correct
            }

            SAMRecordComparator recordComparator = sortOrder.getComparatorInstance();
            if (recordComparator == null) {
                // Should not happen if sort order is coordinate or queryname
                System.err.println("Could not get a comparator for sort order: " + sortOrder + ". Assuming not verifiably sorted.");
                return false;
            }

            SAMRecord previousRecord = null;
            long recordNumber = 0;

            try (CloseableIterator<SAMRecord> iterator = reader.iterator()) {
                while (iterator.hasNext()) {
                    SAMRecord currentRecord = iterator.next();
                    recordNumber++;

                    if (previousRecord != null) {

                        if (recordComparator.fileOrderCompare(previousRecord, currentRecord) > 0) {

                            System.err.println("File is NOT sorted according to " + sortOrder + ".");
                            System.err.println("Violation found at record number " + recordNumber + ":");
                            System.err.println("Previous: " + previousRecord.getReadName() +
                                    (previousRecord.getReadUnmappedFlag() ? " (Unmapped)" :
                                            " Ref: " + previousRecord.getReferenceName() +
                                                    " Start: " + previousRecord.getAlignmentStart()));
                            System.err.println("Current:  " + currentRecord.getReadName() +
                                    (currentRecord.getReadUnmappedFlag() ? " (Unmapped)" :
                                            " Ref: " + currentRecord.getReferenceName() +
                                                    " Start: " + currentRecord.getAlignmentStart()));

                            return false;
                        }
                    }
                    previousRecord = currentRecord;
                }
            }

            //file is empty if there are no record numbers
            if (recordNumber == 0) {
                System.out.println("File is empty. Considered sorted.");
                return true;
            }

            System.out.println("Successfully iterated through " + recordNumber + " records. File appears to be sorted according to " + sortOrder + ".");
            return true;

        }
    }

    public static void convertSamToBam(File samFile, File bamFile, boolean presorted, int compressionLevel) {

        long startTime = System.nanoTime(); //Used to benchmark this method

        if (!samFile.exists()) {
            System.err.println("Error: Input SAM file not found at " + samFile.getAbsolutePath());
            return;
        }

        System.out.println("Starting SAM to BAM conversion (multithreaded compression)...");
        System.out.println("Input SAM: " + samFile.getAbsolutePath());
        System.out.println("Output BAM: " + bamFile.getAbsolutePath());



        try (SamReader samReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(samFile)) {

            SAMFileHeader samHeader = samReader.getFileHeader();

            boolean isPresortedEffective = presorted;
            if (samHeader.getSortOrder() == SAMFileHeader.SortOrder.unsorted) {
                isPresortedEffective = true;
            }

            SAMFileWriterFactory factory = new SAMFileWriterFactory()
                    .setCreateIndex(false)
                    .setCreateMd5File(false)
                    .setCompressionLevel(compressionLevel)
                    .setUseAsyncIo(true);


            try (SAMFileWriter bamWriter = factory.makeBAMWriter(samHeader, isPresortedEffective, bamFile)) {

                long recordCount = 0;
                for (final SAMRecord samRecord : samReader) {
                    bamWriter.addAlignment(samRecord);
                    recordCount++;
                    if (recordCount % 1000000 == 0) {
                        System.out.println("Processed " + recordCount + " records...");
                    }
                }
                long endTime = System.nanoTime();
                System.out.println("Conversion complete. Total records processed: " + recordCount);
                System.out.println("Execution time: " + endTime);
            }

        } catch (IOException e) {
            System.err.println("Error during SAM to BAM conversion: " + e.getMessage());
            e.printStackTrace();
        }


    }

}
