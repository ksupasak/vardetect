/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core;

/**
 *
 * @author soup
 */
import java.util.HashMap;
import java.util.concurrent.LinkedBlockingQueue;

public class ThreadPool {
    private final int nThreads;
    private final PoolWorker[] threads;
    private final LinkedBlockingQueue queue;
    private ThreadPoolCallBack callback;
    private HashMap map;
    private volatile boolean allJobDone;

    public ThreadPool(int nThreads, ThreadPoolCallBack cb) {
        this.nThreads = nThreads;
        this.callback = cb;
        queue = new LinkedBlockingQueue();
        threads = new PoolWorker[nThreads];
        map = new HashMap();
        for (int i = 0; i < nThreads; i++) {
            threads[i] = new PoolWorker();
            threads[i].start();
        }
        allJobDone = false;
    }
    
    

    public void execute(Runnable task, int id) {
        synchronized (queue) {
            map.put(task.hashCode(), id);
            queue.add(task);
            queue.notify();
        }
    }

    void shutdown() {
        System.out.println("Shutdown");
        for (int i = 0; i < nThreads; i++) {
            threads[i] = new PoolWorker();
            threads[i].interrupt();
            
        }
        
        allJobDone = true;
//        System.exit(0);
    }

    public boolean isAllJobDone() {
        return allJobDone;
    }
    
    private class PoolWorker extends Thread {
        public void run() {
            Runnable task;
            
            while (true) {
                synchronized (queue) {
                    while (queue.isEmpty()) {
                        try {
                            queue.wait();
                        } catch (InterruptedException e) {
                            System.out.println("An error occurred while queue is waiting: " + e.getMessage());
                        }
                    }
                    task = (Runnable)queue.poll();
                }
                
                // If we don't catch RuntimeException,
                // the pool could leak threads
                try {
                    task.run();
                } catch (Exception e) {
                    System.err.println("Thread pool is interrupted due to an issue: " + e.getMessage());
                    System.err.println(e.getStackTrace().toString());
                    int id = (int)map.get(task.hashCode());
                    if(callback!=null)callback.finishBatch(id);
                }
                
            }
        }
    }
}
