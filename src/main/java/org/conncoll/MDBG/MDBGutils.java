package org.conncoll.MDBG;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnel;
import com.google.common.hash.Funnels;
import org.jspecify.annotations.Nullable;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.*;


public class MDBGutils {
    // --- SETUP ---
    public static long  expectedKmerCount = 100_000_000; //an estimate of unique k-mers
    public static double  falsePositiveRate = 0.01;      //can tolerate a 1% error rate
    private static final long FNV_OFFSET_BASIS = 0xcbf29ce484222325L;
    private static final long FNV_PRIME        = 0x100000001b3L;
    public record EncodedSequence(byte[] sequence, int validBaseCount) {}


    //compute the 64-bit FNV-1a hash of an ASCII k-mer
    public static long fnv1a64(long kmer) {
        long hash = FNV_OFFSET_BASIS;

        // Process the 8 bytes of the long, from most to least significant
        for (int i = 7; i >= 0; i--) {
            byte b = (byte) (kmer >> (i * 8));
            hash ^= b;
            hash *= FNV_PRIME;
        }
        return hash;
    }

    public static BloomFilter multiThreadedBF(File inputFile, int threads, @Nullable BloomFilter seenOnce) throws IOException {
        long startTime = System.nanoTime(); // for benchmarking time

         BloomFilter<Long> mainBloomFilter = BloomFilter.create(
                Funnels.longFunnel(),
                expectedKmerCount,
                falsePositiveRate
        );


        Queue<BloomFilter<Long>> subBloomFilters = new ConcurrentLinkedQueue<>();

        ExecutorService executorService = Executors.newFixedThreadPool(threads);
        CountDownLatch latch = new CountDownLatch(threads);
        BlockingQueue<String> queue = new LinkedBlockingQueue<>(1000);

        Runnable producerTask = () -> {
            try(BufferedReader reader = new BufferedReader(new FileReader(inputFile))){
                String line;
                while (reader.readLine() !=null){
                    line = reader.readLine(); // actual sequence
                    reader.readLine(); // + stuff
                    reader.readLine(); // quality scores
                    queue.put(line);
                }
            } catch (InterruptedException | IOException e) {
                throw new RuntimeException(e);
            }finally {
                // send in poison pill
                for (int i = 0; i < threads; i++) {
                    try {
                        System.out.println("Producer done reading file");
                        queue.put("HALT");
                    } catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                }
            }
        };

        executorService.submit(producerTask);

        Runnable task;

        if(seenOnce == null){
             task = () ->{
                BloomFilter<Long> localBF = BloomFilter.create(
                        Funnels.longFunnel(),
                        expectedKmerCount,
                        falsePositiveRate
                );

                try {
                    while (true){
                        String sequence = queue.take();
                        if(sequence.equals("HALT")){ //Poison pill
                            break;
                        }
                        // Here is the work each thread will do
                        EncodedSequence encodedSequence = encode(sequence);
                        processKmers(encodedSequence.sequence,encodedSequence.validBaseCount, localBF);
                        //System.out.println("Processed sequence");
                    }
                } catch (InterruptedException e) {
                    System.out.println(Thread.currentThread().getName() + "Finished");
                }finally {
                    subBloomFilters.add(localBF);
                    latch.countDown();
                }
            };

        }else {

             task = () ->{
                BloomFilter<Long> localBF = BloomFilter.create(
                        Funnels.longFunnel(),
                        expectedKmerCount,
                        falsePositiveRate
                );

                try {
                    while (true){
                        String sequence = queue.take();
                        if(sequence.equals("HALT")){ //Poison pill
                            break;
                        }
                        //work thread does
                        EncodedSequence encodedSequence = encode(sequence);
                        processKmersPass2(encodedSequence.sequence,encodedSequence.validBaseCount, localBF, seenOnce);
                        //System.out.println("Processed sequence");
                    }
                } catch (InterruptedException e) {
                    System.out.println(Thread.currentThread().getName() + "Finished");
                }finally {
                    subBloomFilters.add(localBF);
                    latch.countDown();
                }
            };


        }



        for (int i = 0; i < threads; i++) { //for every thread, create a task
            executorService.submit(task);
        }

        try {
            latch.await();//the program will wait here till all threads are done
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }

        executorService.shutdown();
        long endTime = System.nanoTime();

        mergeBloomFilters(mainBloomFilter,subBloomFilters);

        System.out.println("Done with pass 1 and took: " + TimeUnit.NANOSECONDS.toSeconds(endTime-startTime));

        return mainBloomFilter;
    }

    public static EncodedSequence encode(String rawSequence) {
        // A temporary list to hold the 2-bit codes.
        ArrayList<Byte> codes = new ArrayList<>(rawSequence.length());
        for (char rawBase : rawSequence.toCharArray()) {
            switch (rawBase) {
                case 'A': codes.add((byte) 0b00); break;
                case 'C': codes.add((byte) 0b01); break;
                case 'G': codes.add((byte) 0b10); break;
                case 'T': codes.add((byte) 0b11); break;
                // Skips other characters
            }
        }

        int validBaseCount = codes.size();
        int finalArraySize = (validBaseCount + 3) / 4;
        byte[] encodedSequence = new byte[finalArraySize];

        // Now, pack the codes into the single, correctly-sized array.
        for (int i = 0; i < validBaseCount; i++) {
            int byteIndex = i / 4;
            int bitPosition = i % 4;
            int shift = (3 - bitPosition) * 2;
            encodedSequence[byteIndex] |= (codes.get(i) << shift);
        }

        return new EncodedSequence(encodedSequence, validBaseCount);
    }



    public static List<Long> createFastaIndex(File inputFile) throws IOException {
        List<Long> recordStartOffsets = new ArrayList<>();
        long bytesRead = 0;


        try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(inputFile)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {

                    recordStartOffsets.add(bytesRead);
                }
                bytesRead += line.getBytes().length + 1;
            }
        }
        return recordStartOffsets;
    }
    public static List<Long> createFastqIndex(File inputFile) throws IOException {
        List<Long> recordStartOffsets = new ArrayList<>();
        long bytesRead = 0;

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(inputFile)))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("@")) {
                    recordStartOffsets.add(bytesRead);
                }
                bytesRead += line.getBytes().length + 1;
            }
        }
        return recordStartOffsets;
    }

    public static String decodeByteToString(byte packedByte) {
        StringBuilder sb = new StringBuilder(4);

        // An array to map the 2-bit code to a character.
        char[] baseMap = {'A', 'C', 'G', 'T'};

        // Loop four times, once for each base packed in the byte.
        for (int i = 0; i < 4; i++) {
            // The shift determines which 2-bit code we are extracting.
            // 1st base: shift=6, 2nd: shift=4, etc.
            int shift = (3 - i) * 2;

            // Use a mask to isolate the two bits we want.
            // Then, right-shift them to get a value from 0 to 3.
            int code = (packedByte >> shift) & 0b11; // 0b11 is a mask for the last two bits.

            // Append the correct character to our string builder.
            sb.append(baseMap[code]);
        }

        return sb.toString();
    }


    public static void processKmers(byte[] encodedSequence, int validBaseCount, BloomFilter<Long> localBloomFilters) {
        int k = 31;
        if (validBaseCount < k) {
            return; // Sequence is shorter than k.
        }

        // A mask to keep our k-mer at 62 bits (31 bases * 2 bits/base).
        long KMER_MASK = (1L << (k * 2)) - 1;

        // --- Step 1: Manually build the first k-mer ---
        long currentKmer = 0L;
        for (int i = 0; i < k; i++) {
            currentKmer <<= 2; // Make room for the next base.
            byte newBaseCode = getBaseAt(encodedSequence, i); // Get 2-bit code for base i.
            currentKmer |= newBaseCode; // Add it to the k-mer.
        }

        // Process the first k-mer
            localBloomFilters.put(currentKmer);

        // --- Step 2: Slide the window for the rest of the sequence ---
        for (int i = k; i < validBaseCount; i++) {
            // a. Shift the k-mer left by 2 bits, discarding the oldest base.
            currentKmer <<= 2;

            // b. Get the 2-bit code for the new base entering the window.
            byte newBaseCode = getBaseAt(encodedSequence, i);

            // c. Add the new base to the right side of the k-mer.
            currentKmer |= newBaseCode;

            // d. Apply the mask to trim the k-mer back to the correct 62-bit length.
            currentKmer &= KMER_MASK;

            // Process the new overlapping k-mer
            localBloomFilters.put(currentKmer);

        }
    }

    public static void processKmersPass2(byte[] encodedSequence, int validBaseCount, BloomFilter<Long> localBloomFilter, BloomFilter<Long> localSeenOnce) {
        int k = 31;
        if (validBaseCount < k) {
            return; // Sequence is shorter than k.
        }

        // A mask to keep our k-mer at 62 bits (31 bases * 2 bits/base).
        long KMER_MASK = (1L << (k * 2)) - 1;

        // --- Step 1: Manually build the first k-mer ---
        long currentKmer = 0L;
        for (int i = 0; i < k; i++) {
            currentKmer <<= 2; // Make room for the next base.
            byte newBaseCode = getBaseAt(encodedSequence, i); // Get 2-bit code for base i.
            currentKmer |= newBaseCode; // Add it to the k-mer.
        }

        // Process the first k-mer
        if(localSeenOnce.mightContain(currentKmer)){
            localBloomFilter.put(currentKmer);
        }


        // --- Step 2: Slide the window for the rest of the sequence ---
        for (int i = k; i < validBaseCount; i++) {
            // a. Shift the k-mer left by 2 bits, discarding the oldest base.
            currentKmer <<= 2;

            // b. Get the 2-bit code for the new base entering the window.
            byte newBaseCode = getBaseAt(encodedSequence, i);

            // c. Add the new base to the right side of the k-mer.
            currentKmer |= newBaseCode;

            // d. Apply the mask to trim the k-mer back to the correct 62-bit length.
            currentKmer &= KMER_MASK;

            // Process the new overlapping k-mer
            if(localSeenOnce.mightContain(currentKmer)){
                localBloomFilter.put(currentKmer);
            }

        }
    }

    private static void mergeBloomFilters(BloomFilter<Long> outputBloomFilter,Queue<BloomFilter<Long>> subBloomFilters) {

        while (!subBloomFilters.isEmpty()){
            outputBloomFilter.putAll(subBloomFilters.poll());
        }
    }


    public static byte getBaseAt(byte[] encodedSequence, int basePosition) {
        int byteIndex = basePosition / 4;
        int bitPosition = basePosition % 4;
        int shift = (3 - bitPosition) * 2;

        return (byte) ((encodedSequence[byteIndex] >> shift) & 0b11);
    }

    public static void saveBloomFilter(BloomFilter<Long> bf, String path){
        try(FileOutputStream fos = new FileOutputStream(path)){
            bf.writeTo(fos);
            System.out.println("Bloom filter written to: " + path);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static BloomFilter readBloomFilter(String path){
        try(FileInputStream fis = new FileInputStream(path)){
            Funnel<Long> funnel = Funnels.longFunnel();
            BloomFilter bf = BloomFilter.readFrom(fis,funnel);
            return bf;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    public static long reverseComplement(long kmer, int k) {
        long rc = 0;
        for (int i = 0; i < k; i++) {
            // Get the 2 bits for the current base
            long base = kmer & 3;
            // Complement the base (A<->T, C<->G is just XOR with 3)
            // 00(A) ^ 11 = 11(T)
            // 01(C) ^ 11 = 10(G)
            base ^= 3;
            // Shift the reverse complement to make space for the new base
            rc <<= 2;
            // Add the new base
            rc |= base;
            // Move to the next base in the original k-mer
            kmer >>= 2;
        }
        return rc;
    }

    public static String reverseComplement(String sequence) {
        StringBuilder sb = new StringBuilder();
        for (int i = sequence.length() - 1; i >= 0; i--) {
            char base = sequence.charAt(i);
            switch (base) {
                case 'A': sb.append('T'); break;
                case 'T': sb.append('A'); break;
                case 'C': sb.append('G'); break;
                case 'G': sb.append('C'); break;
                default:  sb.append(base); break; // Keep non-ACGT characters as they are
            }
        }
        return sb.toString();
    }


    public static long getCanonical(long kmer, int k) {
        long rc = reverseComplement(kmer, k);
        return kmer < rc ? kmer : rc;
    }

    public static String decode(long kmer, int k) {
        StringBuilder sb = new StringBuilder(k);
        char[] baseMap = {'A', 'C', 'G', 'T'};

        for (int i = 0; i < k; i++) {
            // Get the last 2 bits of the k-mer
            int code = (int) (kmer & 0b11);
            // Append the corresponding base character
            sb.append(baseMap[code]);
            // Shift the k-mer to process the next base
            kmer >>= 2;
        }

        // The loop builds the string in reverse order, so we reverse it.
        return sb.reverse().toString();
    }





}
