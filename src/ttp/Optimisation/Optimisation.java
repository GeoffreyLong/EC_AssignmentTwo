package ttp.Optimisation;


import java.io.*;
import java.util.Comparator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ttp.TTPInstance;
import ttp.TTPSolution;
import ttp.Utils.DeepCopy;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author wagner
 */
public class Optimisation {
    

    public static TTPSolution simpleHeuristic(TTPInstance instance, int[] tour, int maxRuntime) {
    	double[] D = new double[instance.numberOfNodes];
    	double dSum = 0;
    	D[instance.numberOfNodes-1] = 0; 
    	for (int i = instance.numberOfNodes-2; i >= 0; i--) {
    		dSum += instance.distances(tour[i+1], tour[i]);
    		D[tour[i]] = dSum;
    	}
    	double noItemTime = dSum/instance.maxSpeed;
    	double v = (instance.maxSpeed - instance.minSpeed) / instance.capacityOfKnapsack;
    	
    	double[] score = new double[instance.numberOfItems];
    	double[] threshScore = new double[instance.numberOfItems];
    	for (int i = 0; i < instance.numberOfItems; i++) {
    		int cityIdx = instance.items[i][3];
    		double itemCarryTime = D[cityIdx] / (instance.maxSpeed - v*instance.items[i][2]);
    		double itemCycleTime = noItemTime - D[cityIdx] + itemCarryTime;
    		score[i] = instance.items[i][1] - instance.rentingRatio * itemCarryTime;
    		threshScore[i] = instance.rentingRatio*noItemTime 
    				+ (instance.items[i][1] - instance.rentingRatio * itemCycleTime);
    	}
    	
    	// Form array of item indexes to sort
    	Integer[] itemIdx = new Integer[instance.numberOfItems];
    	for (int i = 0; i < itemIdx.length; i++) {
    		itemIdx[i] = instance.items[i][0];
    	}
    	
    	
    	// Heuristic sort
    	Arrays.sort(itemIdx, new Comparator<Integer>() {
    		@Override
    		public int compare(Integer o1, Integer o2) {
    			double diff = score[o1] - score[o2];
    			return (diff == 0) ? 0 : ((diff > 0) ? 1 : -1);
    		}
    	});
    	
    	// Construct solution
    	int[] packingPlan = new int[instance.numberOfItems];
    	int Wc = 0;
    	
    	
    	for (int i = 0; i < instance.numberOfItems; i++) {
    		// If we're not full
    		if ( ((Wc + instance.items[itemIdx[i]][2]) < instance.capacityOfKnapsack) && threshScore[itemIdx[i]] > 0) {
    			//packingPlan[itemIdx[i]] = 1;
    			
    			int itemsPCity=(int)Math.round((double)instance.numberOfItems/(instance.numberOfNodes-1));
    			int cityIndex=instance.items[itemIdx[i]][3];
    			int itemNumber=(int)Math.floor((double)(instance.items[itemIdx[i]][0])/(instance.numberOfNodes-1));
    			
    			int ppIndex = ((cityIndex-1)*itemsPCity)+itemNumber;
    			if(cityIndex<5)
    			System.out.println(ppIndex + " CI " + (cityIndex+1) + " id: " + (instance.items[itemIdx[i]][0]));
    			packingPlan[ppIndex] = 1;
    			
    			Wc += instance.items[itemIdx[i]][2];
    			//System.out.printf("ItemID: %d\n", itemIdx[i]);
    		}
    		if (Wc == instance.capacityOfKnapsack) {
    			break;
    		}
    			//System.out.println(Wc);
    	}
    	System.out.println("Our wc: " +Wc);
    	
    	
   
    	System.out.println(Arrays.toString(packingPlan));
    	TTPSolution s = new TTPSolution(tour, packingPlan);
    	instance.evaluate(s);	
        return s;
    }
    
    public static TTPSolution hillClimber(TTPInstance instance, int[] tour, 
            int mode, 
            int durationWithoutImprovement, int maxRuntime) {
        ttp.Utils.Utils.startTiming();
        
        TTPSolution s = null;
        boolean debugPrint = !true;
        
        int[] packingPlan = new int[instance.numberOfItems];
        
                
        boolean improvement = true;
        double bestObjective = Double.NEGATIVE_INFINITY;
        
        long startingTimeForRuntimeLimit = System.currentTimeMillis()-200;
        
        int i = 0;
        int counter = 0;
        while(counter<durationWithoutImprovement) {
            
            if (i%10==0 /*do the time check just every 10 iterations, as it is time consuming*/
                    && (System.currentTimeMillis()-startingTimeForRuntimeLimit)>=maxRuntime)
                break;
            
            
            if (debugPrint) {
                System.out.println(" i="+i+"("+counter+") bestObjective="+bestObjective); 
            }
            int[] newPackingPlan = (int[])DeepCopy.copy(packingPlan);
            
            boolean flippedToZero = false;
            
            switch (mode) {
                case 1: 
                    // flip one bit
                    int position = (int)(Math.random()*newPackingPlan.length);
//                    newPackingPlan[position] = Math.abs(newPackingPlan[position]-1);
                    if (newPackingPlan[position] == 1) {
                                newPackingPlan[position] = 0;
                                // investigation: was at least one item flipped to zero during an improvement?
//                                flippedToZero = true;
                    } else {
                        newPackingPlan[position] = 1;
                    }
                    break;
                case 2:
                    // flip with probability 1/n
                    for (int j=0; j<packingPlan.length; j++) {
                        if (Math.random()<1d/packingPlan.length)
                            if (newPackingPlan[j] == 1) {
                                newPackingPlan[j] = 0;
                                // investigation: was at least one item flipped to zero during an improvement?
//                                flippedToZero = true;
                            } else {
                                newPackingPlan[j] = 1;
                            }
                    }
                    break;
            }
            
            
            
//            ttp.Utils.Utils.startTiming();
            TTPSolution newSolution = new TTPSolution(tour, newPackingPlan);
            instance.evaluate(newSolution);
//            System.out.println(ttp.Utils.Utils.stopTiming());
            
                        
            /* replacement condition:
             *   objective value has to be at least as good AND
             *   the knapsack cannot be overloaded
             */
            if (newSolution.ob >= bestObjective && newSolution.wend >=0 ) {
                
                // for the stopping criterion: check if there was an actual improvement 
                if (newSolution.ob > bestObjective && newSolution.wend >=0) {
                    improvement = true;
                    counter = 0;
                }
                
                packingPlan = newPackingPlan;
                s = newSolution;
                bestObjective = newSolution.ob;
                
            } else {
                improvement = false;
                counter ++;
            }
            
            i++;
            
        }
        
        long duration = ttp.Utils.Utils.stopTiming();
        s.computationTime = duration;
        return s;
    }
    
    public static int[] linkernTour(TTPInstance instance) {
        int[] result = new int[instance.numberOfNodes+1];
        
        boolean debugPrint = !true;

        String temp = instance.file.getPath();
        int index = temp.indexOf("_");
        String tspfilename = temp;
        if (index==-1) index = tspfilename.indexOf(".");
        String tspresultfilename = System.getProperty("user.dir") + "/" + temp.substring(0,index)+".linkern.tour";
        if (debugPrint) System.out.println("LINKERN: "+tspfilename);
    
        File tspresultfile = new File(tspresultfilename);
        
        
        try {
            if (!tspresultfile.exists()) {
                List<String> command = new ArrayList<String>();
                command.add("./linkern");
                command.add("-o");
                command.add(tspresultfilename);
                command.add(tspfilename);
//                printListOfStrings(command);
                
                ProcessBuilder builder = new ProcessBuilder(command);
                builder.redirectErrorStream(true);
                final Process process = builder.start();
                InputStream is = process.getInputStream();
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);
                String line;
                while ((line = br.readLine()) != null) {
                    if (debugPrint) System.out.println("<LINKERN> "+line);
                }
                
                if (debugPrint) System.out.println("Program terminated?");    
                int rc = process.waitFor();
                if (debugPrint) System.out.println("Program terminated!");
            }
            
            BufferedReader br = new BufferedReader( new FileReader(tspresultfilename));
            // discard the first line
            br.readLine();
            String line = null; 
            for (int i=0; i<result.length-1; i++) {
                line = br.readLine();
                if (debugPrint) System.out.println("<TOUR> "+line);
                index = line.indexOf(" ");
                int number = Integer.parseInt(line.split("\\s+")[0]);
                result[i] = number; 
            }
            
            if (debugPrint) System.out.println(Arrays.toString(result));
          
            
            } catch (Exception ex) {
            	ex.printStackTrace();
            }
        return result;
    }
    
    public static void doAllLinkernTours() {
        
        boolean debugPrint = false;
        
        File f = new File("instances/tsplibCEIL");
        try {
            if (debugPrint) System.out.println(f.getCanonicalPath());
        } catch (IOException ex) {
        }
        
        File[] fa = f.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                boolean result = false;
//                if (name.contains(".ttp") 
                if (name.contains(".tsp") 
                        ) result = true;
                return result;
            }});
        
        if (debugPrint)
            for (File temp:fa) {
                System.out.println(temp.getAbsolutePath());
            }
        
        // create a nonsense instance just to be able to run linkernTour/1 on it
//        TTPInstance instance = new TTPInstance(new File("."));        
//        int[] tour = new int[0];
//        tour = Optimisation.linkernTour(instance);
        
        
        
//        int[] result = new int[instance.numberOfNodes+1];
//        
//        boolean debugPrint = !true;
//        
//        String temp = instance.file.getAbsolutePath();
//        int index = temp.indexOf("_");
        for(File tsp:fa) {
            String tspfilename = tsp.getAbsolutePath();
            int index = tspfilename.indexOf("_");
            if (index==-1) index = tspfilename.indexOf(".");
            String tspresultfilename = tspfilename.substring(0, index) +".linkern.tour";
//            int index = tspfilename.indexOf(".tsp");
//            String tspresultfilename = tspfilename.substring(0, index) +".linkern.tour";
//            String tspresultfilename = tspfilename+".linkern.tour";

            if (debugPrint) System.out.println("LINKERN: "+tspfilename);

            File tspresultfile = new File(tspresultfilename);

            try {
                if (! tspresultfile.exists()) {
                    List<String> command = new ArrayList<String>();
                    command.add("./linkern");
                    command.add("-o");
                    command.add(tspresultfilename);
                    command.add(tspfilename);
//                    printListOfStrings(command);

                    ProcessBuilder builder = new ProcessBuilder(command);
                    builder.redirectErrorStream(true);
                    
                    ttp.Utils.Utils.startTiming();
                    
                    final Process process = builder.start();
                    InputStream is = process.getInputStream();
                    InputStreamReader isr = new InputStreamReader(is);
                    BufferedReader br = new BufferedReader(isr);
                    String line;
                    while ((line = br.readLine()) != null) {
                        if (debugPrint) System.out.println("<LINKERN> "+line);
                    }
                    if (debugPrint) System.out.println("Program terminated?");    
                    int rc = process.waitFor();
                    
                    long duration = ttp.Utils.Utils.stopTiming();
                    
                    System.out.println( new File(tspresultfilename).getName() +" "+duration);
                    
                    if (debugPrint) System.out.println("Program terminated!");
                }
                
                
                
                
                } catch (Exception ex) {
                }
        }
        
    }
    
    public static void printListOfStrings(List<String> list) {
        String result = "";
        for (String s:list)
            result+=s+" ";
        System.out.println(result);
    }
}
