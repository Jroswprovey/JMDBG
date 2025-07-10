package org.conncoll;

import java.util.concurrent.TimeUnit;

import static org.conncoll.Menu.*;

public class ProgressDisplay implements Runnable{
    public volatile boolean isFinished = false;

    @Override
    public void run() {
        long startTime = System.nanoTime();
        long peakMem = 0;

        //Progress text stuff
        while (!isFinished){

            long used = (long) ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) * 0.000001);
            double ratio = (double) totalMinimizers / totalKmers * 100;
            long curTime = TimeUnit.SECONDS.convert(System.nanoTime() - startTime, TimeUnit.NANOSECONDS);

            if (peakMem< used){
                peakMem = used;
            }
            System.out.printf(
                    "\rMemory usage: %dMB │ Peak memory: %dMB  │  Minimizers: %d out of %d  │  Selected: %.2f%%  │  Minimizers in mem test %d  │ Elapsed: %ds ",
                    used,peakMem, totalMinimizers, totalKmers, ratio, minimizerToId.size(), curTime
            );

            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
