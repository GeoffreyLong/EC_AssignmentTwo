package ttp.Optimisation;


import ga.Crossover;
import ga.Mutation;
import ga.Population;
import ga.Config;
import ga.Selection;
import ga.Selection.SelectionType;

import java.awt.Point;
import java.awt.geom.Point2D;
import java.io.*;
import java.util.Collections;
import java.util.Comparator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.javaml.core.kdtree.KDTree;
import ttp.TTPInstance;
import ttp.TTPSolution;
import ttp.Utils.DeepCopy;
import ttp.newrep.City;
import ttp.newrep.Individual;
import ttp.newrep.Item;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author wagner
 */
public class Optimisation {
    
	public static TTPSolution cosolver(TTPInstance instance, int[] tour, int maxRuntime,int H) {
		ttp.Utils.Utils.startTiming();
		Config config = Config.getInstance();
		config.setTtpInstance(instance);
		double[] d = new double[instance.numberOfNodes];
		double[] W = new double[instance.numberOfNodes];
		double[] W2 = new double[instance.numberOfNodes];
		int[] tourRet = new int[tour.length];
		int[] tourDash = new int[tour.length];
		W[0]=0;
		W2[0]=0;
		int[] packingPlanRet = new int[instance.numberOfItems];
		int[] packingPlanDash = new int[instance.numberOfItems];
		double P = Double.NEGATIVE_INFINITY;
		double PDash = Double.NEGATIVE_INFINITY;
		long runtime = 0;
		long startTime = System.currentTimeMillis();
		Individual individual = instance.createIndividual(tour);
		
		//initial Weightings from input tour
		for(int i=0;i < individual.tour.length; i++){
			W[i+1]=individual.tour[i].getWeight();//by tour order
			W2[individual.tour[i].cityId]=individual.tour[i].getWeight();//by city index order
		}
		
		//initial distances from input tour
		d[0]=instance.distances(0, individual.tour[0].cityId);
		for(int i = 0; i < individual.tour.length-1; i++){
			d[i+1] = instance.distances(individual.tour[i].cityId,individual.tour[i+1].cityId);
		}
		d[instance.numberOfNodes-1]=instance.distances(individual.tour[individual.tour.length-1].cityId,0);
		
		while (runtime<maxRuntime){			
			
			packingPlanDash = solveKRP(instance,d,W,individual,H);
			individual=instance.createIndividual(tour, packingPlanDash);
			
			for(int i=0;i < individual.tour.length; i++){
				W[i+1]=individual.tour[i].getWeight();//by tour order
				W2[individual.tour[i].cityId]=individual.tour[i].getWeight();//by city index order
			}
			
			tourDash = solveTSKP(W2,individual);

			PDash = instance.evaluate(individual);
			
			if (PDash>P){
				P=PDash;
				tourRet=tourDash;
				packingPlanRet=packingPlanDash;
				
				d[0]=instance.distances(0, individual.tour[0].cityId);
				for(int i = 0; i < individual.tour.length-1; i++){
					d[i+1] = instance.distances(individual.tour[i].cityId,individual.tour[i+1].cityId);
				}
				d[instance.numberOfNodes-1]=instance.distances(individual.tour[individual.tour.length-1].cityId,0);
			} else {
				break;
			}
			runtime=System.currentTimeMillis()-startTime;
		}
				
		TTPSolution s = new TTPSolution(tourRet, packingPlanRet);
		long duration = ttp.Utils.Utils.stopTiming();
	    s.computationTime = duration;
    	instance.evaluate(s);
        return s;
        
    }
	

	private static int[] solveKRP(TTPInstance instance, double[] d, double[] W, Individual individual,int H){
		
		int[] packingPlanRet = new int[instance.numberOfItems];

		double Pdash = Double.NEGATIVE_INFINITY;
		double P = Double.NEGATIVE_INFINITY;
		double profit = 0;
		double t = 0;	
		int itemsPerCity = instance.numberOfItems / individual.tour.length;
		
		long maxRuntime=100;
		long runtime = 0;
		long startTime = System.currentTimeMillis();
		
		//build initial PP solution
		//could use SH
		TTPSolution s = ppGreedyRegardTour(instance, instance.getTour(individual));
		packingPlanRet=s.packingPlan;
		s.altPrint();
		Individual individualNew=instance.createIndividual(instance.getTour(individual), packingPlanRet);
		
		while (runtime<maxRuntime){	
			// modify the packing plan towards optimal
			//calc profit
			System.out.println("VALUES MUST BE EQUAL: "+individualNew.tour[0].items.size()+" : "+itemsPerCity);
			for (int i = 0; i < individualNew.tour.length; i++){
				for(int j = 0; j < individualNew.tour[i].items.size(); j++){
					//if(packingPlan[(i*itemsPerCity + j)]==1)
					if(individualNew.tour[i].items.get(j).isSelected){
						profit += individualNew.tour[i].items.get(j).profit;
					}
				}
			}
			
			//calc renting rate
			for (int i = 0; i < instance.numberOfNodes; i++){
				t += d[i]/(instance.maxSpeed - W[i]*((instance.maxSpeed-instance.minSpeed)/instance.capacityOfKnapsack));
			}
			
			Pdash = profit - instance.rentingRatio*t;
			
			if (Pdash>P){
				packingPlanRet=instance.getPackingPlan(individualNew);
				P=Pdash;
			}
			
			runtime=System.currentTimeMillis()-startTime;
		}
		
		return packingPlanRet;
	}
	private static int[] solveTSKP(double[] W, ttp.newrep.Individual ind) {
		Config config = Config.getInstance();
		config.setTSKPw(W);
		config.setGenerationMix(true);
		config.setParentSelectionType(SelectionType.ELITISM);
		config.setCrossoverChance(1);
		config.setMutationChance(1);
		int populationSize = 50;
		config.setPopulationSize(populationSize);
		config.setInverOverProbability(0.02);
		config.setTournamentSize(2);
		// set inverOver probability and fitness function
		Mutation mutation = new Mutation(new double[]{0,0,0,0,0, 1});
		Crossover crossover = new Crossover(new double[]{1,0,0,0});
		Selection selection = new Selection(SelectionType.ELITISM);
		Population population = new Population(populationSize-48, ind.tour.length);
		
		ga.Individual currentSol = new ga.Individual();
		currentSol.genotype = new ArrayList<Object>(ind.tour.length);
		for (int i = 0; i < ind.tour.length; i++) {
			currentSol.genotype.add(Integer.toString(ind.tour[i].cityId));
		}
		
		for (int i = 0; i < 48; i++){
			population.population.add(currentSol.clone());
		}
		//population.population.add(currentSol);
		
		
		int numberOfGenerations = 0;
		int maxGeneration = 1000;
		double bestSolution = Double.NEGATIVE_INFINITY;
		ga.Individual bestSolInd=null;
		
		System.out.println("--------------------------------------------------------------------------");
		System.out.println("GEN #     ITER BEST (  POP BEST,    POP AVG,  POP WORST), TIME SINCE ITER START ( OVERALL TIME AVG, OVERALL TIME SUM)");
		while(numberOfGenerations < maxGeneration | true) {
			Population offspring = population.clone();
			
			//offspring = crossover.cross(offspring);
			population = mutation.inverOver(offspring);
			/*
			offspring = mutation.mutate(offspring);
			
			if (config.generationMix){
				population.population.addAll(offspring.population);
			}
			else{
				population.population = offspring.population;
			}
			population = selection.select(population);
			*/
			numberOfGenerations++;
			
			/// calc data store best worst and avg
			//bestF,avgF,worstF,bestInd,worstInd
			Double[] data = population.getStats();

			if (1/data[0] > bestSolution){
				bestSolution = 1/data[0];
				bestSolInd=population.population.get(data[3].intValue());
			}
			
			System.out.println("G: "+String.format("%5d",numberOfGenerations)+" "+String.format("%10.3f",bestSolution) + " ("+String.format("%10.2f",1/data[0])+", "+String.format("%10.2f",1/data[1])+", "+String.format("%10.2f",1/data[2])+"), \n");
		}
		
		int[] sol = new int[bestSolInd.genotype.size()+2];
		sol[0] = 0;
		for (int i = 0; i < bestSolInd.genotype.size(); i++) {
			sol[i+1] = Integer.parseInt((String)bestSolInd.genotype.get(i));
		}
		sol[sol.length-1] = 0;
		
		return sol;
}
	
    public static TTPSolution simpleHeuristic(TTPInstance instance, int[] tour, int maxRuntime) {
    	ttp.Utils.Utils.startTiming();
    	double[] D = new double[instance.numberOfNodes];
    	double dSum = 0;
    	//D[instance.numberOfNodes-1] = 0; 
    	D[0] = 0;
    	//for (int i = instance.numberOfNodes-2; i >= 0; i--) {
    	for (int i = tour.length-2; i >= 0; i--) { // if >= D[0] is set to total distance
    		dSum += instance.distances(tour[i+1], tour[i]);
    		D[tour[i]] = dSum;
    	}
    	double noItemTime = dSum/instance.maxSpeed;
    	double v = (instance.maxSpeed - instance.minSpeed) / instance.capacityOfKnapsack;
    	
    	final double[] score = new double[instance.numberOfItems];
    	double[] threshScore = new double[instance.numberOfItems];
    	for (int i = 0; i < instance.numberOfItems; i++) {
    		int cityIdx = instance.items[i][3];
    		double itemCarryTime = D[cityIdx] / (instance.maxSpeed - v*instance.items[i][2]);
    		double itemCycleTime = noItemTime - D[cityIdx] + itemCarryTime;
    		score[instance.items[i][0]] = instance.items[i][1] - instance.rentingRatio * itemCarryTime;
    		threshScore[instance.items[i][0]] = instance.rentingRatio*noItemTime 
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
    			return (diff == 0) ? 0 : ((diff > 0) ? -1 : 1); // swapped order
    		}
    	});
    	
    	// Construct solution
    	int[] packingPlan = new int[instance.numberOfItems];
    	int Wc = 0;
    	
    	TTPSolution s = new TTPSolution(tour, packingPlan);
    	instance.evaluate(s);	
    	for (int i = 0; i < instance.numberOfItems; i++) {
    		// If we're not full
    		//if ( ((Wc + instance.items[itemIdx[i]][2]) < instance.capacityOfKnapsack) && threshScore[itemIdx[i]] > 0) {
    		if ((Wc + instance.items[itemIdx[i]][2]) < instance.capacityOfKnapsack) {
    			int arrIndex=-1;
    			int itemsPCity=(int)Math.round((double)instance.numberOfItems/(instance.numberOfNodes-1));
    			int cityIndex=instance.items[itemIdx[i]][3];
    			int itemNumber=(int)Math.floor((double)(instance.items[itemIdx[i]][0])/(instance.numberOfNodes-1));
    			
    			for (int j = 1; j<tour.length; j++){
    				if (tour[j]==cityIndex){
    					arrIndex=j-1;    					
    					break;
    				}	
    			}

    			int ppIndex = (arrIndex*itemsPCity)+itemNumber;
    			
    			// Only use the item if it produces a better solution
    			packingPlan[ppIndex] = 1;
    			TTPSolution newSol = new TTPSolution(tour,packingPlan);
    			instance.evaluate(newSol);
    			if (newSol.ob > s.ob) {
    				s = newSol;
        			Wc += instance.items[itemIdx[i]][2];
    			} else {
    				packingPlan[ppIndex] = 0;
    			}
    			
    			//System.out.println("i: "+arrIndex+" ppI: "+ppIndex + " CI " + (cityIndex) + " id: " + (instance.items[itemIdx[i]][0])+" IN: "+itemNumber+" WGC: "+ instance.items[itemIdx[i]][2] +" WGT: "+Wc);
    			//System.out.printf("ItemID: %d\n", itemIdx[i]);
    		}
    		if (Wc == instance.capacityOfKnapsack) {
    			break;
    		}
    	}
    	//System.out.println("Our wc: " +Wc);
    	
    	
   
    	//System.out.println(Arrays.toString(packingPlan));
        long duration = ttp.Utils.Utils.stopTiming();
        s.computationTime = duration;
        return s;
    }
    
    /**
     * Sorts the items by profit over cost ratio
     * It then selects items greedily (highest profit cost ratio first)
     * It changes the packing plan for the corresponding object as long as the fitness increases
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseTwoSolutionTwo(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime, boolean debug){
    	ttp.Utils.Utils.startTiming();
    	long startingTimeForRunLimit = System.currentTimeMillis();
    	
    	List<double[]> items = getProfitWeightRatios(tour, instance, startingTimeForRunLimit, maxRuntime);

        TTPSolution solution = exerciseTwoSolutionTwoLogic(instance, tour, durationWithoutImprovement, maxRuntime, items, startingTimeForRunLimit, debug);
        
        solution.computationTime = ttp.Utils.Utils.stopTiming();
        
    	return solution;
    }
    
    /**
     * Same as the other solutionTwo, except instead of profit cost ratio it is by the cutoff weight
     * It then selects items greedily (highest profit cost ratio first)
     * It changes the packing plan for the corresponding object as long as the fitness increases
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseTwoSolutionTwoAlt(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime, boolean debug){
    	ttp.Utils.Utils.startTiming();
    	long startingTimeForRunLimit = System.currentTimeMillis();
                
    	List<double[]> items = getWeightCutoffs(tour, instance, startingTimeForRunLimit, maxRuntime);
    	TTPSolution solution = exerciseTwoSolutionTwoLogic(instance, tour, durationWithoutImprovement, maxRuntime, items, startingTimeForRunLimit, debug);
         
    	solution.computationTime = ttp.Utils.Utils.stopTiming();
    	
    	return solution;
    }
    
    /**
     * Same as previous two, but instead of costly sorting simply iterate through items in reverse tour order
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @param debug
     * @return
     */
    public static TTPSolution exerciseTwoSolutionTwoAltTwo(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime, boolean debug){
    	ttp.Utils.Utils.startTiming();
    	long startingTimeForRunLimit = System.currentTimeMillis();
                
    	List<double[]> items = new LinkedList<double[]>();
    	int itemsPerCity = instance.numberOfItems / (tour.length - 2);
    	for (int i = tour.length-1; i >= 0; i--){
    		
    		int cityIndex = tour[i];
    		// Do not want cityIndex of 0, 
    		// this node has no items and will cause array out of bounds on itemIndex lookup
    		if (cityIndex != 0){
				for (int j = 0; j < itemsPerCity; j++){
					int itemIndex = (tour.length-2) * j + cityIndex-1;
					items.add(new double[]{itemIndex,0});
				}
    		}
    	}
    	
    	TTPSolution solution = exerciseTwoSolutionTwoLogic(instance, tour, durationWithoutImprovement, maxRuntime, items, startingTimeForRunLimit, debug);
         
    	solution.computationTime = ttp.Utils.Utils.stopTiming();
    	
    	return solution;
    }
    
    /**
     * No sorting at all, no linked list structure
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @param debug
     * @return
     */
    public static TTPSolution exerciseTwoSolutionTwoAltThree(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime, boolean debug){
    	long startingTimeForRunLimit = System.currentTimeMillis();
    	int[] packingPlan = new int[instance.numberOfItems];
		int[] packingPlanClone = packingPlan.clone();
		TTPSolution newSolution = new TTPSolution(tour, packingPlan);
	    instance.evaluate(newSolution);
	    double lastSolutionOb = -Double.MAX_VALUE;
	    int didNotImprove = 0;
	    
	    int k = 0;
	    
	    // Break the loop when the individual has not improved within a certain number of iterations
	    while(didNotImprove <= durationWithoutImprovement){
	    	// Counter for improvement logic
			if (lastSolutionOb == newSolution.ob){
				didNotImprove ++;
			}
			else{
				didNotImprove = 0;
			}
	
			// Only change the packing plan if the solution is valid
			if (newSolution.wend >= 0){
				packingPlanClone = packingPlan.clone();
			}
			
			if(debug) System.out.println(lastSolutionOb);
			lastSolutionOb = newSolution.ob;
			k++;
			if (k%10==0 /*do the time check just every 10 iterations, as it is time consuming*/
	                && (System.currentTimeMillis()-startingTimeForRunLimit)>=maxRuntime)
	            break;
			
			// Iterate through the items array
			for (int i = 0; i < instance.numberOfItems; i ++){
				
				// Greedily choose the highest ratio not yet encountered
				
				// If this item has not yet been picked up
				// Pick it up and check if this new packing plan improves the fitness
				//		if it does then keep this new packing plan
				// 		else revert the plan back
				if (packingPlan[i] 	== 0){
					packingPlan[i] = 1;
	    			TTPSolution tempSolution = new TTPSolution(tour, packingPlan);
	                instance.evaluate(tempSolution);
	                
	                if (tempSolution.ob > newSolution.ob){
	                	break;
	                }
	                else{
	                	packingPlan[i] = 0;
	                }
				}
				
				// If the item has been picked up
				// Drop the item and check if this new packing plan improves the fitness
				//		if it does then keep this new packing plan
				//		else revert the plan back
				if (packingPlan[i] == 1){
					packingPlan[i] = 0;
	    			TTPSolution tempSolution = new TTPSolution(tour, packingPlan);
	                instance.evaluate(tempSolution);
	                
	                if (tempSolution.ob > newSolution.ob){
	                	break;
	                }
	                else{
	                	packingPlan[i] = 1;
	                }
				}
				
			}
			newSolution = new TTPSolution(tour, packingPlan);
	        instance.evaluate(newSolution);
		}
	    
	    // Create a solution with the packing plan clone
	    // The packing plan clone is guaranteed to be a valid solution (wend >= 0)
	    // Whereas the packing plan may be invalid
	    TTPSolution solution = new TTPSolution(tour, packingPlanClone);
	    instance.evaluate(solution);
	    
		return solution;
	}
    
    public static TTPSolution exerciseTwoSolutionTwoLogic(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime, List<double[]> items, long startingTimeForRunLimit, boolean debug){
    	if(debug) System.out.println("Sorting took " + (System.currentTimeMillis()-startingTimeForRunLimit));
    	int[] packingPlan = new int[instance.numberOfItems];
    	int[] packingPlanClone = packingPlan.clone();
    	TTPSolution newSolution = new TTPSolution(tour, packingPlan);
        instance.evaluate(newSolution);
        double lastSolutionOb = -Double.MAX_VALUE;
        int didNotImprove = 0;
        
        int k = 0;
        
        // Break the loop when the individual has not improved within a certain number of iterations
        while(didNotImprove <= durationWithoutImprovement){
        	// Counter for improvement logic
    		if (lastSolutionOb == newSolution.ob){
    			didNotImprove ++;
    		}
    		else{
    			didNotImprove = 0;
    		}

    		// Only change the packing plan if the solution is valid
    		if (newSolution.wend >= 0){
    			packingPlanClone = packingPlan.clone();
    		}
    		
    		if(debug) System.out.println(lastSolutionOb);
    		lastSolutionOb = newSolution.ob;
    		k++;
    		if (k%10==0 /*do the time check just every 10 iterations, as it is time consuming*/
                    && (System.currentTimeMillis()-startingTimeForRunLimit)>=maxRuntime)
                break;
    		
    		// Iterate through the items array
    		for (int i = 0; i < instance.numberOfItems; i ++){
    			
    			// Greedily choose the highest ratio not yet encountered
    			int index = (int) items.get(i)[0];
    			
    			// If this item has not yet been picked up
    			// Pick it up and check if this new packing plan improves the fitness
    			//		if it does then keep this new packing plan
    			// 		else revert the plan back
    			if (packingPlan[index] 	== 0){
    				packingPlan[index] = 1;
	    			TTPSolution tempSolution = new TTPSolution(tour, packingPlan);
	                instance.evaluate(tempSolution);
	                
	                if (tempSolution.ob > newSolution.ob){
	                	break;
	                }
	                else{
	                	packingPlan[index] = 0;
	                }
    			}
    			
    			// If the item has been picked up
    			// Drop the item and check if this new packing plan improves the fitness
    			//		if it does then keep this new packing plan
    			//		else revert the plan back
    			if (packingPlan[index] == 1){
    				packingPlan[index] = 0;
	    			TTPSolution tempSolution = new TTPSolution(tour, packingPlan);
	                instance.evaluate(tempSolution);
	                
	                if (tempSolution.ob > newSolution.ob){
	                	break;
	                }
	                else{
	                	packingPlan[index] = 1;
	                }
    			}
    			
    		}
    		newSolution = new TTPSolution(tour, packingPlan);
            instance.evaluate(newSolution);
    	}
        
        // Create a solution with the packing plan clone
        // The packing plan clone is guaranteed to be a valid solution (wend >= 0)
        // Whereas the packing plan may be invalid
        TTPSolution solution = new TTPSolution(tour, packingPlanClone);
        instance.evaluate(solution);
        
		return solution;
    }
    
    public static TTPSolution randomLinkernTours(TTPInstance instance, long maxRunTime) {
    	ttp.Utils.Utils.startTiming();
    	long startTime = System.currentTimeMillis();
    	int numIterations = 0;
    	long elapsedTime = 0;
    	TTPSolution bestSol = null;
    	double bestObj = Double.NEGATIVE_INFINITY;
    	boolean two = false;
    	boolean five = false;
    	boolean ten = false;
    	do {
    		// Generate a new linkern tour
        	int[] linTour = linkernTour(instance.file.getPath(), instance.numberOfNodes+1);
    		//int[] linTour = linkernTour(instance);
        	// Generate a packing plan
        	int itemsPerCity = instance.numberOfItems / (instance.numberOfNodes - 1);
        	TTPSolution sol = flipTourCheck(instance, linTour);
        	instance.evaluate(sol);    	    
        	
    		elapsedTime = System.currentTimeMillis() - startTime;
    		if (elapsedTime >= 60000*2 && !two) {
    			System.out.println("Two Min: " + bestObj);
    			two = true;
    		}
    		if (elapsedTime >= 60000*5 && !five) {
    			System.out.println("Five Min: " + bestObj);
    			five = true;
    		}
    		if (elapsedTime >= 60000*10 && !ten) {
    			System.out.println("Ten Min: " + bestObj);
    			ten = true;
    		}
        	
        	if (sol.ob > bestObj) {
        		bestSol = sol;
        		bestObj = sol.ob;
        	}
    		//System.out.printf("BestObj: %f, NewObj: %f\n",bestObj, sol.ob);
    		numIterations++;
    	} while (elapsedTime <= maxRunTime);
    	bestSol.computationTime = elapsedTime;
    	return bestSol;
    }
    
    public static TTPSolution exerciseThreeLinkernCycles(TTPInstance instance) {
    	ttp.Utils.Utils.startTiming();
    	// Find best items in terms of profit/weight
        double[] profitWeightRatio = new double[instance.numberOfItems];

		for(int i = 0; i < instance.numberOfItems; i++){
			profitWeightRatio[i] = instance.items[i][1] / instance.items[i][2];
		}
		
		// Sort items
		double[][] sortData = new double[instance.numberOfItems][2];
		
		for(int i = 0; i<instance.numberOfItems; i++){
			sortData[i][0]=i;
			sortData[i][1]=profitWeightRatio[i];
		}
		
		Comparator<double[]> newComp = new Comparator<double[]>() {
    		@Override
    		public int compare(double[] s1, double[] s2) {
    			return -Double.compare(s1[1], s2[1]);
		    }
    	};
    	Arrays.sort(sortData,newComp);
    	
    	
    	// Find best city indexes
    	int maxCities = instance.numberOfNodes/2; // The item solution must contain at most this number of cities
    	List<Integer> bestCities = new ArrayList<Integer>();
    	Set<Integer> bestCitiesSet = new HashSet<Integer>();
    	double totalWeight = 0;
    	int index = 0;
    	while (totalWeight < instance.capacityOfKnapsack && index < instance.numberOfItems) {
    		int sortedItemIndex = (int)sortData[index][0];
    		if ((totalWeight + instance.items[sortedItemIndex][2]) <= instance.capacityOfKnapsack) {
    			if (!bestCitiesSet.contains(instance.items[sortedItemIndex][3])) {
    				bestCities.add(instance.items[sortedItemIndex][3]);
    				bestCitiesSet.add(instance.items[sortedItemIndex][3]);
    				if (bestCitiesSet.size() >= maxCities) {
    					break;
    				}
    			}
    			totalWeight += instance.items[sortedItemIndex][2];
    		}
    		index++;
    	}
    	System.out.println("BEST CITIES: " + bestCities);
    	// Generate linkern compatible tsp file with best cities 
    	try {
    		BufferedWriter writer = new BufferedWriter(new FileWriter("linkernCycles.txt"));
    		String output = "";
    		output += "PROBLEM NAME: linkernCycles\n";
    		output += "DIMENSION: " + (bestCities.size()+1) + "\n";
    		output += "EDGE_WEIGHT_TYPE: CEIL_2D\n";
    		output += "NODE_COORD_SECTION	(INDEX, X, Y): \n";
    		output += String.format("%d \t %d %d\n", 0, (int)instance.nodes[0][1], (int)instance.nodes[0][2]);
    		for (int cityId : bestCities) {
    			output += String.format("%d \t %d %d\n", cityId, (int)instance.nodes[cityId][1], (int)instance.nodes[cityId][2]);
    		}
    		writer.write(output);
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    		return null;
    	}
    	
    	// Generate linkern between best items from start city
    	int[] bestCitiesTour = linkernTour("linkernCycles.txt", bestCities.size()+2);
    	// Tour indicies converted to city order in file, map back to normal range
    	for (int i = 1; i < bestCitiesTour.length-1; i++) {
    		bestCitiesTour[i] = bestCities.get( bestCitiesTour[i] - 1);
    	}
    	System.out.println("Linkern for item cities: " + Arrays.toString(bestCitiesTour));
    	
    	// 'Fill' in the generated tour
    	// Construct kd tree for nearest neighbour lookups
    	KDTree tree = new KDTree(2);
    	for (int i = 0; i < instance.numberOfNodes; i++) {
    		tree.insert(new double[]{instance.nodes[i][1], instance.nodes[i][2]}, i);
    	}
    	//System.out.println( tree.nearest(new double[]{288,149}));
    	List<Integer> filledTour = new ArrayList<Integer>();
    	Set<Integer> filledTourSet = new HashSet<Integer>();
    	filledTour.add(0);
    	filledTourSet.add(0);
    	int numNearestNeighbours = 40;
    	int prevCity = 0;
    	for (int i = 1; i < bestCitiesTour.length; i++) {
    		int currentCity =  bestCitiesTour[i];
    		double remainingDistance = instance.distances(prevCity,currentCity);
    		// fill in path between prev city and current city, only using nodes that get us closer to currentCity
    		while (prevCity != currentCity) {
	    		Object[] nn = tree.nearest(new double[]{instance.nodes[prevCity][1],instance.nodes[prevCity][2]},numNearestNeighbours);
	    		boolean found = false;
	    		// get line equation: Ax + By + C = 0
	    		double x1 = instance.nodes[prevCity][1];
	    		double x2 = instance.nodes[currentCity][1];
	    		double y1 = instance.nodes[prevCity][2];
	    		double y2 = instance.nodes[currentCity][2];
	    		double A = y1 - y2;
	    		double B = x2 - x1;
	    		double C = x1*y2 - x2*y1;
	    		

	    		// perp dist threshold
	    		double perpDistThreshold = Math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	    		perpDistThreshold /= 0.5;
	    		
	    		for (Object o : nn) {
	    			double newRemainingDistance = instance.distances((int)o,currentCity);
	    			double perpDist = Math.abs(A*instance.nodes[(int)o][1] + B*instance.nodes[(int)o][2] + C);
	    			perpDist /= Math.sqrt(A*A + B*B);
	    			// && perpDist < perpDistThreshold
	    			if ((newRemainingDistance < remainingDistance && !filledTourSet.contains((int)o) && !bestCitiesSet.contains((int)o)) || ((int)o) == currentCity) {
	    				prevCity = (int)o;
	    				filledTour.add(prevCity);
	    				filledTourSet.add(prevCity);
		    			remainingDistance = newRemainingDistance;
		    			found = true;
	    				break;
	    			}
	    		}
	    		if (!found) {
	    			filledTour.add(currentCity);
	    			filledTourSet.add(currentCity);
	    			prevCity = currentCity;
	    		}
    		}
    	}
    	
    	
    	// Generate linkern tour from the set of cities not including the above
    	List<Integer> remainingCities = new ArrayList<Integer>();
    	for (int i = 0; i < instance.numberOfNodes; i++) {
    		if (!filledTourSet.contains(i)) {
    			remainingCities.add(i);
    		}
    	}
    	System.out.println("REMAINING CITIES: " + remainingCities);
    	// Generate linkern compatible tsp file with best cities 
    	try {
    		BufferedWriter writer = new BufferedWriter(new FileWriter("linkernCycles.txt"));
    		String output = "";
    		output += "PROBLEM NAME: linkernCycles\n";
    		output += "DIMENSION: " + (remainingCities.size()+1) + "\n";
    		output += "EDGE_WEIGHT_TYPE: CEIL_2D\n";
    		output += "NODE_COORD_SECTION	(INDEX, X, Y): \n";
    		output += String.format("%d \t %d %d\n", 0, (int)instance.nodes[0][1], (int)instance.nodes[0][2]);
    		for (int cityId : remainingCities) {
    			output += String.format("%d \t %d %d\n", cityId, (int)instance.nodes[cityId][1], (int)instance.nodes[cityId][2]);
    		}
    		writer.write(output);
    		writer.close();
    	} catch (IOException e) {
    		e.printStackTrace();
    		return null;
    	}
    	
    	// Generate linkern for remaining cities
    	int[] remainingCitiesTour = linkernTour("linkernCycles.txt", remainingCities.size()+2);
    	// Tour indices converted to city order in file, map back to normal range
    	for (int i = 1; i < remainingCitiesTour.length-1; i++) {
    		remainingCitiesTour[i] = remainingCities.get( remainingCitiesTour[i] - 1);
    	}
    	System.out.println("Linkern for remaining cities: " + Arrays.toString(remainingCitiesTour));
    	
    	// Combine linkern cycles with best cities tour at end
    	int[] combinedTour = new int[instance.numberOfNodes+1];
    	for (int i = 0; i < remainingCitiesTour.length; i++) {
    		combinedTour[i] = remainingCitiesTour[i];
    	}
    	for (int i = remainingCitiesTour.length - 1, j = 0; i < combinedTour.length; i++, j++) {
    		combinedTour[i] = filledTour.get(j+1);
    	}
    	System.out.println("COMBINED TOUR: " + Arrays.toString(combinedTour));
    	
    	TTPSolution sol = ppGreedyRegardTour(instance, combinedTour);
    	instance.evaluate(sol);    	    	
        long duration = ttp.Utils.Utils.stopTiming();
        sol.computationTime = duration;
        
        try {
        	BufferedWriter bw = new BufferedWriter(new FileWriter("filledCitiesTour"));
        	for (Integer i : filledTour) {
        		bw.write(String.format("%d\t%d\t%d\n",i, (int)instance.nodes[i][1], (int)instance.nodes[i][2]));
        	}
        	bw.close();
        	bw = new BufferedWriter(new FileWriter("bestCitiesTour"));
        	for (Integer i : bestCitiesTour) {
        		bw.write(String.format("%d\t%d\t%d\n",i, (int)instance.nodes[i][1], (int)instance.nodes[i][2]));
        	}
        	bw.close();
        	bw = new BufferedWriter(new FileWriter("combinedTour"));
        	for (Integer i : combinedTour) {
        		bw.write(String.format("%d\t%d\t%d\n",i, (int)instance.nodes[i][1], (int)instance.nodes[i][2]));
        	}
        	bw.close();
        	bw = new BufferedWriter(new FileWriter("remainingTour"));
        	for (Integer i : remainingCitiesTour) {
        		bw.write(String.format("%d\t%d\t%d\n",i, (int)instance.nodes[i][1], (int)instance.nodes[i][2]));
        	}
        	bw.close();
        } catch (Exception e) {
        	e.printStackTrace();
        }
        
        
    	return sol;
    }
        
    
    
    public static List<double[]> getProfitWeightRatios(int[] tour, TTPInstance instance, long startingTimeForRunLimit, int maxRuntime){
    	double [] distances = new double[tour.length];
		double tourDistance = 0;
		
		// Get the distance from the node to the end
		// Create array of cityId -> distance to end of tour for easy lookup
		java.awt.geom.Point2D.Double lastPoint = new Point.Double(instance.nodes[tour[0]][1], instance.nodes[tour[0]][2]);
		// Iterate through tour backwards
		for (int i = tour.length-1; i>=0; i--){
			int index = tour[i];
			java.awt.geom.Point2D.Double point = new Point.Double(instance.nodes[index][1], instance.nodes[index][2]);
    		tourDistance += point.distance(lastPoint);
    		
    		distances[index] = tourDistance;
		}    	
    		
		// List will store the index and profit/cost ratio of all items
		List<double[]> items = new LinkedList<double[]>();

		int itemsPerCity = instance.numberOfItems / (tour.length - 2);
		int runCount = 0;
		for (int i = tour.length-1; i >= 0; i--){
    		runCount++;
    		if (runCount%10==0 /*do the time check just every 10 iterations, as it is time consuming*/
                    && (System.currentTimeMillis()-startingTimeForRunLimit)>=maxRuntime)
                break;
    		
    		int cityIndex = tour[i];
    		// Do not want cityIndex of 0, 
    		// this node has no items and will cause array out of bounds on itemIndex lookup
    		if (cityIndex != 0){
				for (int j = 0; j < itemsPerCity; j++){
					int itemIndex = (tour.length-2) * j + cityIndex-1;
	    			int[] item = instance.items[itemIndex];
	    			int itemWeight = item[2];
	    			int itemProfit = item[1];
	    			
	    			double cost = instance.rentingRatio * distances[cityIndex];
	    			
	    			// Takes the symbolic place of the other weights that make up the total weight in the knapsack
	    			// Perhaps make it ratio of avg weight : avg distance?
	    			// TODO This is temp penalty, make a better one
	    			double penalty = distances[cityIndex] + itemWeight / 2;
	    			cost /= instance.maxSpeed - (itemWeight + penalty) * (instance.maxSpeed - instance.minSpeed) / instance.capacityOfKnapsack;
	    			
	    			double ratio = itemProfit / cost;
	    			
	    			// Save the itemIndex and the ratio into a double array
	    			// Both of these pieces of data are necessary in the sort
	    			double[] nodeArray = new double[2];
	    			nodeArray[0] = (i-1)*itemsPerCity + j;
		    		nodeArray[1] = ratio;
		    		
		    		// Add item to the list according to its ratio (descending)
		    		for (int k = 0; k <= items.size(); k++){
		    			if (k == items.size()){
		    				items.add(nodeArray);
		    				break;
		    			}
		    			else{
			    			if (nodeArray[1] >= items.get(k)[1]){
			    				items.add(k, nodeArray);
			    				break;
			    			}
		    			}
		    		}
				}
			}
    	}
    	return items;
    }
    
    public static List<double[]> getWeightCutoffs(int[] tour, TTPInstance instance, long startingTimeForRunLimit, int maxRuntime){
    	double [] distances = new double[tour.length];
		double tourDistance = 0;
		
		// Get the distance from the node to the end
		// Create array of cityId -> distance to end of tour for easy lookup
		java.awt.geom.Point2D.Double lastPoint = new Point.Double(instance.nodes[tour[0]][1], instance.nodes[tour[0]][2]);
		// Iterate through tour backwards
		for (int i = tour.length-1; i>=0; i--){
			int index = tour[i];
			java.awt.geom.Point2D.Double point = new Point.Double(instance.nodes[index][1], instance.nodes[index][2]);
    		tourDistance += point.distance(lastPoint);
    		
    		distances[index] = tourDistance;
		} 
		
		int itemsPerCity = instance.numberOfItems / (tour.length - 2);
		List<double[]> items = new LinkedList<double[]>();
		int runCount = 0;
    	for (int i = 0; i <= tour.length - 1; i++){
    		runCount++;
    		if (runCount%10==0 /*do the time check just every 10 iterations, as it is time consuming*/
                    && (System.currentTimeMillis()-startingTimeForRunLimit)>=maxRuntime)
                break;
    		int cityIndex = tour[i];
    		// Do not want cityIndex of 0, 
    		// this node has no items and will cause array out of bounds on itemIndex lookup
    		if (cityIndex != 0){
				for (int j = 0; j < itemsPerCity; j++){
					int itemIndex = (tour.length-2) * j + cityIndex-1;
	    			int[] item = instance.items[itemIndex];
	    			int itemWeight = item[2];
	    			int itemProfit = item[1];
	    			
	    			
	    			// Remember to edit this, need a wt on the cost without picking up the item
	    			double cut = (instance.maxSpeed - instance.rentingRatio * distances[cityIndex] / (itemProfit + instance.rentingRatio*distances[cityIndex]/instance.maxSpeed)) * instance.capacityOfKnapsack;
	    			cut /= (instance.maxSpeed - instance.minSpeed);
	    			cut -= itemWeight;
	    			
	    			// Save the itemIndex and the ratio into a double array
	    			// Both of these pieces of data are necessary in the sort
	    			double[] nodeArray = new double[3];
	    			nodeArray[0] = (i-1)*itemsPerCity + j;
		    		nodeArray[1] = cut;
		    		nodeArray[2] = itemWeight;
		    		
		    		/*
		    		// Add item to the list according to its ratio (descending)
		    		for (int k = 0; k <= items.size(); k++){
		    			if (k == items.size()){
		    				items.add(nodeArray);
		    				break;
		    			}
		    			else{
			    			if (nodeArray[1] >= items.get(k)[1]){
			    				items.add(k, nodeArray);
			    				break;
			    			}
		    			}
		    		}
		    		*/
		    		
		    		// Ascending positive followed by ascending negative
		    		if (nodeArray[1] <= 0){
		    			for (int k = 0; k <= items.size(); k++){
			    			if (k == items.size()){
			    				items.add(nodeArray);
			    				break;
			    			}
			    			else{
				    			if (nodeArray[1] >= items.get(k)[1]){
				    				items.add(k, nodeArray);
				    				break;
				    			}
			    			}
			    		}
		    		}
		    		else{
		    			for (int k = 0; k <= items.size(); k++){
			    			if (k == items.size()){
			    				items.add(nodeArray);
			    				break;
			    			}
			    			else{
				    			if (nodeArray[1] <= items.get(k)[1] || items.get(k)[1] <= 0){
				    				items.add(k, nodeArray);
				    				break;
				    			}
			    			}
			    		}
		    		}
				}
    		}
    	}
		return items;
    }
    
    
    
    
    /**
     * If the tour is fixed, this should output the cutoff weight for  the knapsack
     * If the knapsack is already this weight, theoretically adding the item will not add profit
     * 
     * @param instance
     * @param tour
     * @return
     */
    public static TTPSolution getWeightCuts(TTPInstance instance, int[] tour){
    	double [] distances = new double[tour.length];
		double tourDistance = 0;
		
		// Get the distance from the node to the end
		// Create array of cityId -> distance to end of tour for easy lookup
		java.awt.geom.Point2D.Double lastPoint = new Point.Double(instance.nodes[tour[0]][1], instance.nodes[tour[0]][2]);
		// Iterate through tour backwards
		for (int i = tour.length-1; i>=0; i--){
			int index = tour[i];
			java.awt.geom.Point2D.Double point = new Point.Double(instance.nodes[index][1], instance.nodes[index][2]);
    		tourDistance += point.distance(lastPoint);
    		
    		distances[index] = tourDistance;
		} 
		
		int itemsPerCity = instance.numberOfItems / (tour.length - 2);
		List<double[]> items = new LinkedList<double[]>();
    	for (int i = 0; i <= tour.length - 1; i++){
    		int cityIndex = tour[i];
    		// Do not want cityIndex of 0, 
    		// this node has no items and will cause array out of bounds on itemIndex lookup
    		if (cityIndex != 0){
				for (int j = 0; j < itemsPerCity; j++){
					int itemIndex = (tour.length-2) * j + cityIndex-1;
	    			int[] item = instance.items[itemIndex];
	    			int itemWeight = item[2];
	    			int itemProfit = item[1];
	    			
	    			
	    			// Remember to edit this, need a wt on the cost without picking up the item
	    			double cut = (instance.maxSpeed - instance.rentingRatio * distances[cityIndex] / (itemProfit + instance.rentingRatio*distances[cityIndex]/instance.maxSpeed)) * instance.capacityOfKnapsack;
	    			cut /= (instance.maxSpeed - instance.minSpeed);
	    			cut -= itemWeight;
	    			
	    			// Save the itemIndex and the ratio into a double array
	    			// Both of these pieces of data are necessary in the sort
	    			double[] nodeArray = new double[3];
	    			nodeArray[0] = i*itemsPerCity + j - 1;
		    		nodeArray[1] = cut;
		    		nodeArray[2] = itemWeight;
		    		items.add(nodeArray);
		    		
		    		// Add item to the list according to its ratio (descending)
		    		for (int k = 0; k <= items.size(); k++){
		    			if (k == items.size()){
		    				items.add(nodeArray);
		    				break;
		    			}
		    			else{
			    			if (nodeArray[1] >= items.get(k)[1]){
			    				items.add(k, nodeArray);
			    				break;
			    			}
		    			}
		    		}
				}
    		}
    	}

    	for (int k = 0; k < items.size(); k++){ 		
			System.out.println(items.get(k)[1] + " " + items.get(k)[0] + " " + items.get(k)[2]);
		}

    	return null;
    }

     /**
     * Might actually be the same as Two New, but it uses a resetter for the packing plan... 

     * @param instance
     * @param tour
     * @param durationWithoutImprovement -> Run with a value of around five
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseThreeSolutionTwo(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
    	ttp.Utils.Utils.startTiming();
    	long startingTimeForRunLimit = System.currentTimeMillis();
    	Individual individual = instance.createIndividual(tour);
        int didNotImprove = 0;
        TTPSolution solution = null;
        int[] packingPlan = new int[instance.numberOfItems];
        double oldFitness = 0;
        int tempCount = 0;
        double bestFitness = 0;
        
        int[] bestPlan = packingPlan.clone();
        int[] bestTour = tour.clone();
		int k=0;
        while(didNotImprove <= durationWithoutImprovement){
    		int randNodeIndex = (int) (Math.random() * individual.tour.length);
    		
			k++;
			if (k%100==0 /*do the time check just every 10 iterations, as it is time consuming*/
	                && (System.currentTimeMillis()-startingTimeForRunLimit)>=maxRuntime)
	            break;
    		
    		solution = new TTPSolution(tour, packingPlan);
    		instance.evaluate(solution);
    		if (oldFitness == solution.ob){
    			tempCount ++;
    		}
    		else {
    			tempCount = 0;
    			oldFitness = solution.ob;
    		}
    		if (oldFitness > bestFitness){
    			bestPlan = packingPlan.clone();
    			bestTour = tour.clone();
        		//System.out.println(bestFitness);
        		didNotImprove = 0;
    			bestFitness = oldFitness;
    		}
    		
    		if (didNotImprove == durationWithoutImprovement){
    			// Note that it still usually improves beyond this point, however, it takes a while to get here
    			// Would probably be best to run this on a timer instead
    			break;
    		}
    		
			if (randNodeIndex > 1 && randNodeIndex < individual.tour.length){
	    		for (int i = 0; i < 5; i++){
	    			int[] tourNew = tour.clone();
	    			int[] packingPlanNew = packingPlan.clone();
	    			
	    			Individual tempInd = instance.createIndividual(tourNew, packingPlanNew);
	    			double randSwap = Math.random();
	    			if (randSwap <= 0.1){
	    				City cityTemp = tempInd.tour[randNodeIndex - 1];
	    				tempInd.tour[randNodeIndex - 1] = tempInd.tour[randNodeIndex];
	    				tempInd.tour[randNodeIndex] = cityTemp;
	    			}
	    			
	    			double packingSwap = Math.random();
	    			for (Item item : tempInd.tour[randNodeIndex - 1].items){
		    			if (packingSwap <= 0.1){
		    				if (item.isSelected) item.isSelected = false;
		    				else item.isSelected = true;
		    			}
	    			}
	    			for (Item item : tempInd.tour[randNodeIndex].items){
		    			if (packingSwap <= 0.1){
		    				if (item.isSelected) item.isSelected = false;
		    				else item.isSelected = true;
		    			}
	    			}
	    			
	    			TTPSolution newSolution = new TTPSolution(instance.getTour(tempInd), instance.getPackingPlan(tempInd));
	    	        instance.evaluate(newSolution);
	    	        // 1000 is a random parameter, but it seems to work well
	    	        if (tempCount >= 1000){
	    	        	System.out.println("revamp");
	    	        	packingPlan = new int[packingPlan.length];
	    	        	tempCount = 0;
	    	        	didNotImprove ++;
	    	        	break;
	    	        }
	    			if (newSolution.ob > solution.ob && newSolution.wend >= 0){
	    				tour = instance.getTour(tempInd).clone();
	    	    		packingPlan = instance.getPackingPlan(tempInd).clone();
	    			}
	    		}
			}
        }
        solution = new TTPSolution(bestTour, bestPlan);
        instance.evaluate(solution);
        solution.computationTime = ttp.Utils.Utils.stopTiming();
    	return solution;
    } 
        
    public static TTPSolution flipTourCheck(TTPInstance instance, int[] tour) {
        int[] tourOrig= new int[tour.length];
        tourOrig=tour.clone();
        TTPSolution newSolution = null;
        
    	for(int i = 1; i<tour.length/2;i++){ //perform the flip of a tour
        	int temp = tour[i];
        	tour[i]=tour[tour.length-1-i];
        	tour[tour.length-1-i]=temp;
    	}
        
    	//apply packing plan and apply result
    	TTPSolution origSolution = Optimisation.ppGreedyRegardTour(instance, tourOrig);
    	instance.evaluate(origSolution);
    	TTPSolution flipSolution = Optimisation.ppGreedyRegardTour(instance, tour);
    	instance.evaluate(flipSolution);
    	if(origSolution.ob>flipSolution.ob){
    		newSolution=origSolution;
    	}else{
    		newSolution=flipSolution;
    	}

        instance.evaluate(newSolution);
        return newSolution;
    }
    
    /**
     * Greedy Heuristic based packing plan builder, with no regard to tour
     * 
     * @param instance
     * @param tour
     * @param individual
     * @param H
     * @return TTPSolution
     */
    public static TTPSolution ppGreedyDisregardTour(TTPInstance instance, int[] tour) {
        ttp.Utils.Utils.startTiming();

        int[] packingPlan = new int[instance.numberOfItems];
        int[] cityIndex = new int[instance.numberOfNodes];
        int[] cityTourIndex = new int[instance.numberOfItems];
        
        double[] weights = new double[instance.numberOfItems];
        double[] weightsRatio = new double[instance.numberOfItems];
        double totW = 0;
        double[] profits = new double[instance.numberOfItems];
        double[] profitsRatio = new double[instance.numberOfItems];
        double totP = 0;
        double[] values = new double[instance.numberOfItems];
        double[] dSoFar = new double[instance.numberOfNodes];
        double[] dToGo = new double[instance.numberOfNodes];
        
        double MAXWEIGHT = instance.capacityOfKnapsack;
        
        int itemsPerCity = instance.numberOfItems / (instance.numberOfNodes-1);
		dSoFar[0]=0;
		cityTourIndex[0]=0;

		for(int i = 1; i < tour.length-1; i++){
			dSoFar[i] = dSoFar[i-1]+instance.distances(tour[i-1],tour[i]);
			cityIndex[tour[i]]=i;
		}
        
		for(int i = 0; i < instance.numberOfNodes; i++){
			dToGo[i] = dSoFar[instance.numberOfNodes-1]-dSoFar[i];
		}
		
		for(int i = 0; i < instance.numberOfItems; i++){
			cityTourIndex[i]=cityIndex[instance.items[i][3]];
			totW+=instance.items[i][2];
			weights[i]=instance.items[i][2];
			totP+=instance.items[i][2];
			profits[i]=instance.items[i][1];
		}

		
		for(int i = 0; i < instance.numberOfItems; i++){
			weightsRatio[i]=weights[i]/totW;
			profitsRatio[i]=profits[i]/totP;
			values[i]=profitsRatio[i]/weightsRatio[i];
		}
		

		//add the items to the PP
		double totalWeight = 0;
		
		double[][] sortData = new double[instance.numberOfItems][2];
		
		for(int i = 0; i<instance.numberOfItems; i++){
			sortData[i][0]=i;
			sortData[i][1]=values[i];
		}
		
		Comparator<double[]> newComp = new Comparator<double[]>() {
    		@Override
    		public int compare(double[] s1, double[] s2) {
    			return -Double.compare(s1[1], s2[1]);
		    }
    	};
    	
    	//System.out.println("Sorting "+instance.numberOfItems+" items...");
    	Arrays.sort(sortData,newComp);

    	//System.out.println("Filling Packing Plan");
		int index=0;
		
		Set<Integer> visitedCities = new HashSet<Integer>();
		double jump=Math.ceil(instance.numberOfItems/20);
		while(totalWeight<MAXWEIGHT && index<instance.numberOfItems && jump>=2){
			int bestValueIndex=(int)sortData[index][0];
			//System.out.println(instance.items[bestValueIndex][3]);
			//add it as long as it doesn't break capacity			
			if(totalWeight+weights[bestValueIndex]<=instance.capacityOfKnapsack){
				int ppIndex=(cityTourIndex[bestValueIndex]-1)*itemsPerCity + (int)(bestValueIndex/(tour.length-2));
				packingPlan[ppIndex]=1;
				totalWeight+=weights[bestValueIndex];
				visitedCities.add(tour[cityTourIndex[bestValueIndex]]);
			}
			index++;			
		}

        long duration = ttp.Utils.Utils.stopTiming();
        //System.out.println("TIME TAKEN: "+duration+" .. Total Weight: "+totalWeight);
        TTPSolution newSolution = new TTPSolution(tour, packingPlan);
        newSolution.computationTime = duration;
        instance.evaluate(newSolution);
        return newSolution;
    }
    
    /**
     * Greedy Heuristic based packing plan builder, with regard to tour
     * 
     * @param instance
     * @param tour
     * @param individual
     * @param H
     * @param jump 
     * @return TTPSolution
     */
    public static TTPSolution ppGreedyRegardTour(TTPInstance instance, int[] tour) {
        ttp.Utils.Utils.startTiming();

        double jump = instance.numberOfItems/100;
        
        int[] packingPlan = new int[instance.numberOfItems];
        int[] cityIndex = new int[instance.numberOfNodes];
        int[] cityTourIndex = new int[instance.numberOfItems];
        
        double[] profitVSweight = new double[instance.numberOfItems];

        double[] weights = new double[instance.numberOfItems];
        double[] weightsRatio = new double[instance.numberOfItems];
        double totW = 0;
        double[] profits = new double[instance.numberOfItems];
        double[] profitsRatio = new double[instance.numberOfItems];
        double totP = 0;
        double[] values = new double[instance.numberOfItems];
        double[] dSoFar = new double[instance.numberOfNodes];
        double[] dToGo = new double[instance.numberOfNodes];
        
        double MAXWEIGHT = instance.capacityOfKnapsack;
        
        int itemsPerCity = instance.numberOfItems / (instance.numberOfNodes-1);
		dSoFar[0]=0;
		cityTourIndex[0]=0;

		for(int i = 1; i < tour.length-1; i++){
			dSoFar[i] = dSoFar[i-1]+instance.distances(tour[i-1],tour[i]);
			cityIndex[tour[i]]=i;
		}
        
		for(int i = 0; i < instance.numberOfNodes; i++){
			dToGo[i] = dSoFar[instance.numberOfNodes-1]-dSoFar[i];
		}
		
		for(int i = 0; i < instance.numberOfItems; i++){
			cityTourIndex[i]=cityIndex[instance.items[i][3]];
			totW+=instance.items[i][2];
			weights[i]=instance.items[i][2];
			totP+=instance.items[i][2];
			profits[i]=instance.items[i][1];
		}

		
		for(int i = 0; i < instance.numberOfItems; i++){
			weightsRatio[i]=weights[i]/totW;
			profitsRatio[i]=profits[i]/totP;
			profitVSweight[i]=profitsRatio[i]/weightsRatio[i];
		}
		
		for(int i = 0; i < instance.numberOfItems; i++){//this is the important weightings
			if(itemsPerCity==1){
				values[i]=Math.pow(profitVSweight[i],4)+profitsRatio[i]-weightsRatio[i]*instance.rentingRatio*(dToGo[cityTourIndex[i]]);//best for 01s
			}else{
				values[i]=Math.pow(profitVSweight[i]*profits[i],1)/Math.pow(instance.rentingRatio*(dToGo[cityTourIndex[i]]/dToGo[0]),1);//best for 05s
			}	
		}

		//add the items to the PP
		double totalWeight = 0;
		int count=0;

		double lastOB=Double.NEGATIVE_INFINITY;
		
		double[][] sortData = new double[instance.numberOfItems][2];
		
		for(int i = 0; i<instance.numberOfItems; i++){
			sortData[i][0]=i;
			sortData[i][1]=values[i];
		}
		
		Comparator<double[]> newComp = new Comparator<double[]>() {
    		@Override
    		public int compare(double[] s1, double[] s2) {
    			return -Double.compare(s1[1], s2[1]);
		    }
    	};
    	
    	//System.out.println("Sorting "+instance.numberOfItems+" items...");
    	Arrays.sort(sortData,newComp);

    	//System.out.println("Filling Packing Plan");
		int noImprovement=0;
		int index=0;
		
		int[] packingPlanOld = new int[instance.numberOfItems];
		int indexOld = 0;
		double weightOld = 0;
		boolean jumpSet=false;
		if(jump>1){
			jumpSet=true;
		}else{
			jumpSet=false;
			jump=2;
		}
		
		while(totalWeight<MAXWEIGHT && count<instance.numberOfItems && jump>=2){
			//System.out.println(100*(index/instance.numberOfItems)+"%");
			int bestValueIndex=(int)sortData[index][0];
			
			count++;
			//add it as long as it doesn't break capacity			
			if(totalWeight+weights[bestValueIndex]<=instance.capacityOfKnapsack){
				int ppIndex=(cityTourIndex[bestValueIndex]-1)*itemsPerCity + (int)(bestValueIndex/(tour.length-2));
				packingPlan[ppIndex]=1;
				totalWeight+=weights[bestValueIndex];
				//System.out.println("I: "+bestValueIndex+" .. P: "+profits[bestValueIndex]+" .. W: "+weights[bestValueIndex]+" .. C: "+cityTourIndex[bestValueIndex]+"/"+(tour.length-2)+" ... V: "+values[bestValueIndex]+" PvW: "+profitVSweight[bestValueIndex]);
				
				if(!jumpSet || index%jump==0){
					
					TTPSolution s = new TTPSolution(tour, packingPlan);
			        instance.evaluate(s);
			        //System.out.println(jump+" .. "+s.ob+" .<?. "+lastOB+" .. "+index+" .. "+noImprovement);
			        if(s.ob<lastOB){//remove if id doesn't improve OB
			        	totalWeight-=weights[bestValueIndex];
			        	packingPlan[ppIndex]=0;
			        	noImprovement++;
			        	
			        	if(jumpSet){
				        	packingPlan=packingPlanOld.clone();
				        	index=indexOld;
				        	totalWeight=weightOld;
			        		jump=Math.ceil(jump/2);
			        	}

			        	//System.out.println("WORSE: "+jump+" .. "+s.ob+" .. "+index+" .. "+noImprovement);
			        }else{
			        	weightOld=totalWeight;
			        	indexOld=index;
			        	noImprovement=0;
			        	lastOB=s.ob;
			        	packingPlanOld=packingPlan.clone();
			        	//System.out.println("BETTER: "+jump+" .. "+s.ob+" .. "+index+" .. "+noImprovement);
					}	
				}					        
			}
			index++;			
		}
		
        long duration = ttp.Utils.Utils.stopTiming();
        //System.out.println("TIME TAKEN: "+duration+" .. Total Weight: "+totalWeight);
        TTPSolution newSolution = new TTPSolution(tour, packingPlan);
        newSolution.computationTime = duration;
        instance.evaluate(newSolution);
        return newSolution;
    }
    
    public static TTPSolution insertion(TTPInstance instance, int[] tour, int[] packingPlan, int maxRuntime){
    	ttp.Utils.Utils.startTiming();
    	long startTime = System.currentTimeMillis();
    	long elapsedTime=0;
    	Individual individual = instance.createIndividual(tour);
        individual = instance.createIndividual(tour, packingPlan);
        
    	double[] cityProfits = new double[tour.length-2];
    	double[] tourCityProfits = new double[tour.length-2];
    	double[] cityWeights = new double[tour.length-2];
    	double[] tourCityWeights = new double[tour.length-2];
    	double[] tourCityWeightsCumulative = new double[tour.length-2];
    	double[][] itemCosts = new double[instance.numberOfItems][3];
    	
    	double[] distToCity = new double[tour.length-2];

    	for(int i = 0; i<instance.numberOfItems; i++){
    		itemCosts[i][2]=i;
    		itemCosts[i][1]=instance.items[i][3];
    	    itemCosts[i][0]=(double)instance.items[i][1]/instance.items[i][2];
    	}
    	
		Comparator<double[]> newComp = new Comparator<double[]>() {
    		@Override
    		public int compare(double[] s1, double[] s2) {
    			return -Double.compare(s1[0], s2[0]);
		    }
    	};
    	
    	//System.out.println("Sorting "+instance.numberOfItems+" items...");
    	Arrays.sort(itemCosts,newComp);
    	
    	// CALCULATE CITY WEIGHTS AND PROFITS AND DISTANCES TRAVELLED TO CITY
    	for(int i = 0; i<cityWeights.length;i++){
    		cityProfits[individual.tour[i].cityId-1] = individual.tour[i].getEmptyProfit();
    		tourCityProfits[i] = individual.tour[i].getProfit();
    		cityWeights[individual.tour[i].cityId-1] = individual.tour[i].getEmptyWeight();
    		tourCityWeights[i] = individual.tour[i].getWeight();
    		if (i>0){
    			tourCityWeightsCumulative[i] = tourCityWeightsCumulative[i-1]+tourCityWeights[i];
    			distToCity[i]=distToCity[i-1]+instance.distances(i+1, i);
    		}else{
    			tourCityWeightsCumulative[i] = tourCityWeights[i];
    			distToCity[i]=instance.distances(i+1, i);
    		}
    	}
    	
    	double overallDist=distToCity[cityWeights.length-1]+instance.distances(cityWeights.length-1, 0);
    	
    	// WORK OUT WHAT CITYS ARE BEST TO MOVE BASED ON WEIGHTS*DISTANCE TO TRAVEL
    	double[][] movePref = new double[cityWeights.length][2];
    	for(int i = 0; i<cityWeights.length;i++){
    		movePref[i][1]=i;
    		movePref[i][0]=tourCityWeights[i]*(overallDist-distToCity[i]);
    	}
    	
    	Arrays.sort(movePref,newComp);

    	// FOR ALL CITIES THAT ARE SUGGESTED TO BE MOVED TRY AND FIT THEM ELSEWHERE IN THE TOUR, CHECK NEW OB OR HUERISTIC VALUE
    	double moveThresh=1;
    	int betterMoves=0;
    	for(int i = 0; i<cityWeights.length;i++){
    		if(movePref[i][0]<moveThresh){
    			break;
    		}
    		//j=tourIndex
    		int j=(int)movePref[i][1];
    		//System.out.println("Attempting to move city "+individual.tour[j].cityId+" from tour position "+j+". Move Ratio: "+movePref[i][0]);
    		
			double distanceWithoutOld=0;
			double distanceWithOld = 0;
			double removedDistance=0;
    		if(j==individual.tour.length-1){
    			continue;
    		}else if (j==0){
    			distanceWithoutOld =instance.distances(0,individual.tour[j+1].cityId);
    			distanceWithOld = instance.distances(individual.tour[j].cityId,0)+instance.distances(individual.tour[j].cityId,individual.tour[j+1].cityId);    		
    			removedDistance=distanceWithoutOld-distanceWithOld;
    		}else{
    			distanceWithoutOld = instance.distances(individual.tour[j-1].cityId,individual.tour[j+1].cityId);
    			distanceWithOld = instance.distances(individual.tour[j].cityId,individual.tour[j-1].cityId)+instance.distances(individual.tour[j].cityId,individual.tour[j+1].cityId);    		
    			removedDistance=distanceWithoutOld-distanceWithOld;
    		}
    		
    		//System.out.println("Distance W: "+distanceWithOld+" ... D WOUT: "+distanceWithoutOld+" ... RM D:"+removedDistance);
    		
    		// LOOK FOR PLACE TO REFIT CITY    		
    		//only want to add later\
    		double currentValue= Double.POSITIVE_INFINITY;
    		double bestValue= Double.POSITIVE_INFINITY;
    		int bestIndex=-1;
    		
    		//put at ni and push all forward (so basically put after ni)    		
    		for(int ni=j+1;ni<cityWeights.length;ni++){
    			double distanceWithoutNew=0;
    			double distanceWithNew = 0;
    			double addedDistance = 0;
    			if(ni==j){// ni == j can compare with current spot (doesn't move)	
    				addedDistance=removedDistance;
    			}else if(ni==individual.tour.length-1){
        			distanceWithoutNew = instance.distances(individual.tour[ni].cityId,0);
        			distanceWithNew = instance.distances(individual.tour[j].cityId,individual.tour[ni].cityId)+instance.distances(individual.tour[j].cityId,0);    		
        			addedDistance=distanceWithNew-distanceWithoutNew;
        		}else{
        			distanceWithoutNew = instance.distances(individual.tour[ni].cityId,individual.tour[ni+1].cityId);
        			distanceWithNew = instance.distances(individual.tour[j].cityId,individual.tour[ni].cityId)+instance.distances(individual.tour[j].cityId,individual.tour[ni+1].cityId);    		
        			addedDistance=distanceWithNew-distanceWithoutNew;
        		}
    			
    			currentValue=(addedDistance+removedDistance)*(tourCityWeightsCumulative[ni]-tourCityWeights[j])/overallDist;
    			    			
    			if(currentValue<=bestValue){
    				//System.out.println(ni+" .. "+currentValue+" .. "+removedDistance+" .. "+addedDistance);
    				bestValue=currentValue;
    				bestIndex=ni;
    			}
    		}
    		
    		Individual old = instance.createIndividual(instance.getTour(individual), instance.getPackingPlan(individual));
    		TTPSolution oldSolution = new TTPSolution(instance.getTour(old), instance.getPackingPlan(individual));
            instance.evaluate(oldSolution);
            double oldOB = oldSolution.ob;
              		
    		// INSERT CITY LATER DOWN TOUR
    		for(int mi=bestIndex;mi>j;mi--){
    			City temp = individual.tour[mi];
    			individual.tour[mi]=individual.tour[j];
    			individual.tour[j]=temp;
    		}
    		
    		TTPSolution newSolution = new TTPSolution(instance.getTour(individual), instance.getPackingPlan(individual));
            instance.evaluate(newSolution);
            double newOB = newSolution.ob;
            
    		if(newOB>oldOB){
    			System.out.println(j+" ...: "+bestIndex+" .. "+bestValue+" .. "+oldOB+" .. "+newOB);
    			betterMoves+=1;
    		}else{
    			individual=old;
    		}
    		elapsedTime = System.currentTimeMillis() - startTime;
    		if (elapsedTime>maxRuntime){
    			break;
    		}
    		
    	}
    	// IF NEW PLACEMENT IS BETTER, REARRANGE TOUR (DOUBLE CHECK IT HELPS BEFORE OVERWRITING?)
        //System.out.println("TIME TAKEN: "+duration+" .. Total Weight: "+totalWeight);
        TTPSolution newSolution = new TTPSolution(instance.getTour(individual), instance.getPackingPlan(individual));
        newSolution.computationTime = elapsedTime;
        instance.evaluate(newSolution);
        return newSolution;
    }
    
    public static TTPSolution inversion(TTPInstance instance, int[] tour, int[] packingPlan, int maxRuntime){
    	ttp.Utils.Utils.startTiming();
   	
    	Individual individual = instance.createIndividual(tour);
        individual = instance.createIndividual(tour, packingPlan);
        Individual individualOld= instance.createIndividual(tour, packingPlan);
		TTPSolution oldS = new TTPSolution(tour, packingPlan);
        instance.evaluate(oldS);
        
        double[] cityWeights = new double[tour.length];
        double[] fromStart = new double[tour.length];
        double[] fromEnd = new double[tour.length];
        
        for(int i = 1; i<tour.length-1; i++){
        	cityWeights[i]=cityWeights[i-1]+individual.tour[i-1].getWeight();
        	fromStart[i]=fromStart[i-1]+instance.distances(tour[i], tour[i-1]);
        }
        
        for(int i = 1; i<tour.length; i++){
        	fromEnd[i]=fromStart[tour.length-1]-fromStart[i];        
        }
        
        int itemsPerCity = instance.numberOfItems / (instance.numberOfNodes-1);
        
        int betterCount=0;
        int betterCountB=0;
        
    	//for(int dist=tour.length-2;dist>1;dist--){
    	for(int dist=4;dist>1;dist--){
    		
    		for(int a = 1; a+dist<tour.length; a++){
    			int b=a+dist-1;
    			double valueBefore=0;
    			for(int c = a; c<=b; c++){
    				valueBefore+=cityWeights[c]*fromEnd[c];
    			}
    			
    			double valueAfter=0;
    			double[] tempWeight= new double[dist];
    			double[] tempFromEnd= new double[dist];
    			for(int c = a; c<=b; c++){
    				if(c==a){
    					tempWeight[dist-(1+c-a)]=tempWeight[dist-(1+c-a)]+(cityWeights[c]-cityWeights[c-1])+cityWeights[a-1];
    				}else{
    					tempWeight[dist-(1+c-a)]=tempWeight[dist-(1+c-a)]+(cityWeights[c]-cityWeights[c-1]);
    				}
    				tempFromEnd[dist-(1+c-a)]=(fromEnd[b+1]+(fromStart[c]-fromStart[a-1])-instance.distances(tour[a-1], tour[a])+instance.distances(tour[a],tour[b+1]));
    				valueAfter+=tempWeight[dist-(1+c-a)]*tempFromEnd[dist-(1+c-a)];
    			}
    			

    			System.out.println("D: "+dist+" A: "+a+" B: "+b+" BF: "+valueBefore+" AV: "+valueAfter+" BC: "+betterCount+" BCB: "+betterCountB);
    			//if(valueAfter<valueBefore){
    			if(true){
    				
    				//System.out.println("-----------D: "+dist+" A: "+a+" B: "+b+" BF: "+valueBefore+" AV: "+valueAfter+" BC: "+betterCount+" BCB: "+betterCountB);
    				betterCount++;
    				
    				int swaps = (int) Math.floor(dist/2);//how many swap operations
    				//System.out.println(indexA+", "+indexB+", "+swaps);//testing
   		        
    				for (int c = 0; c < swaps; c++) {
    					int temp = tour[a+c];//store temp
    					tour[a+c]=tour[b-c];
    					tour[b-c]=temp;
    				}
   				
    				for (int c = 0; c < swaps; c++) {
    					City t = individual.tour[a+c-1];
    					individual.tour[a+c-1]=individual.tour[b-c-1];
    					individual.tour[b-c-1]=t;
    				}	
    				
    				TTPSolution tempS = new TTPSolution(instance.getTour(individual), instance.getPackingPlan(individual));
    		        instance.evaluate(tempS);
    		        //System.out.println(oldS.ob + " ... N: "+ tempS.ob);
    		        
    		        if(tempS.ob>oldS.ob){
        				for(int c = a; c<=b; c++){
        					cityWeights[c]=tempWeight[c-a];
        					fromEnd[c]=tempFromEnd[c-a];
        				}
    		        	betterCountB++;
    		        	System.out.println(oldS.ob + " ... YY: "+ tempS.ob);
    		        	oldS.ob=tempS.ob;
        		        
    		        	individualOld=instance.createIndividual(instance.getTour(individual),instance.getPackingPlan(individual));
    		        }else{
    		        	individual=instance.createIndividual(instance.getTour(individualOld),instance.getPackingPlan(individualOld));
    		        }
    			}
    		}
    	}
    	
    	// IF NEW PLACEMENT IS BETTER, REARRANGE TOUR (DOUBLE CHECK IT HELPS BEFORE OVERWRITING?)
        long duration = ttp.Utils.Utils.stopTiming();
        //System.out.println("TIME TAKEN: "+duration+" .. Total Weight: "+totalWeight);
        TTPSolution newSolution = new TTPSolution(instance.getTour(individual), instance.getPackingPlan(individual));
        newSolution.computationTime = duration;
        instance.evaluate(newSolution);
        return newSolution;
    }

    /**
     * Currently this is a purely stochastic crossover alg
     * It crosses either the full tour or the full packing plan
     * 
     * Steps
     * 	Select two parents at random
     * 	Create two children nodes
     * 		Mutate
     * 			Swap two nodes right next to each other and invert these nodes' item selection
     * 			The chance of either of these happening is based on a certain probability
     * 		Crossover
     * 			For each pairing (16 of them because pairing of individual and self is possible) of the four individuals
     * 				Create a new child that is the packing plan of one and the tour of another
     * 				If this child is better than the individual occupying one of two of the indices the parents were drawn from
     * 					Replace the parent with the child
     * 
     * NOTE this solution has lots of room for optimization in terms of runtime
     * Also the mutation logic could do with some improvement
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseFourSolutionOne(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
    	int popSize = 10;
    	int noImprove = 0;
    	Individual[] population = new Individual[popSize];
    	int itemsPerCity = instance.numberOfItems/(tour.length - 2);
    	for (int i = 0; i < popSize; i++){
    		population[i] = new Individual(instance.nodes, instance.items, tour, itemsPerCity);
    	}
    	
    	TTPSolution bestSol = null;
    	double bestOb = -Double.MAX_VALUE;
    	double lastOb = 0;
    	
    	while(noImprove < durationWithoutImprovement){
    		if (lastOb == bestOb){
    			noImprove ++;
    		}
    		else{
    			lastOb = bestOb;
    			noImprove = 0;
    		}
    		
    		for (int i = 0; i < popSize; i++){
    			int indIndexOne = (int) (Math.random() * popSize);
    			int indIndexTwo = (int) (Math.random() * popSize);
    			Individual indOne = population[indIndexOne];
    			Individual mutateOne = instance.createIndividual(instance.getTour(indOne), instance.getPackingPlan(indOne));
    			Individual indTwo = population[indIndexTwo];
    			Individual mutateTwo = instance.createIndividual(instance.getTour(indTwo), instance.getPackingPlan(indTwo));
    			Individual[] mutates = {mutateOne, mutateTwo};	
    			Individual[] allInds = {indOne, mutateOne, indTwo, mutateTwo};

        		for (Individual ind : mutates){
        			int randNodeIndex = (int) (Math.random() * indOne.tour.length);
        			double randSwap = Math.random();
        			if (randSwap <= 0.1){
        				City cityTemp = ind.tour[(randNodeIndex - 1 + ind.tour.length) % ind.tour.length];
        				ind.tour[(randNodeIndex - 1 + ind.tour.length) % ind.tour.length] = ind.tour[randNodeIndex];
        				ind.tour[randNodeIndex] = cityTemp;
        			}
        			
        			double packingSwap = Math.random();
	    			for (Item item : ind.tour[(randNodeIndex - 1 + ind.tour.length) % ind.tour.length].items){
		    			if (packingSwap <= 0.1){
		    				if (item.isSelected) item.isSelected = false;
		    				else item.isSelected = true;
		    			}
	    			}
	    			for (Item item : ind.tour[randNodeIndex].items){
		    			if (packingSwap <= 0.1){
		    				if (item.isSelected) item.isSelected = false;
		    				else item.isSelected = true;
		    			}
	    			}
        		}
        		
        		double oddSol = -Double.MAX_VALUE;
        		double evenSol = -Double.MAX_VALUE;
				for(int j = 0; j < allInds.length; j++){
					for(int k = 0; k < allInds.length; k++){
						TTPSolution sol = new TTPSolution(instance.getTour(allInds[j]), instance.getPackingPlan(allInds[k]));
		    			instance.evaluate(sol);
		    			int index = (i+1)*j;
		    			if (sol.wend >= 0){
			    			if (index % 2 == 0 && sol.ob > evenSol){
			    				evenSol = sol.ob;
			    				population[indIndexOne] = instance.createIndividual(instance.getTour(allInds[j]), instance.getPackingPlan(allInds[k]));
			    			}
			    			if (index % 2 == 1 && sol.ob > oddSol){
			    				oddSol = sol.ob;
			    				population[indIndexTwo] = instance.createIndividual(instance.getTour(allInds[j]), instance.getPackingPlan(allInds[k]));
			    			}
			    			if (sol.ob > bestOb){
			    				bestSol = sol;
			    				bestOb = sol.ob;
			    			}
		    			}
					}
        		}
    		}
    	}
    	return bestSol;
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
    
    public static int[] linkernTour(String tourFileName, int numNodes) {
        int[] result = new int[numNodes];
        
        boolean debugPrint = !true;
        File tourFile = new File(tourFileName);
        String temp = tourFile.getPath();
        int index = temp.indexOf("_");
        String tspfilename = temp;
        if (index==-1) index = tspfilename.indexOf(".");
        String tspresultfilename = System.getProperty("user.dir") + "/" + temp.substring(0,index)+".linkern.tour.alt";
        if (debugPrint) System.out.println("LINKERN: "+tspfilename);
  
        
        try {
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
            br.close();
            
            br = new BufferedReader( new FileReader(tspresultfilename));
            // discard the first line
            br.readLine();
            line = null; 
            for (int i=0; i<result.length-1; i++) {
                line = br.readLine();
                if (debugPrint) System.out.println("<TOUR> "+line);
                index = line.indexOf(" ");
                int number = Integer.parseInt(line.split("\\s+")[0]);
                result[i] = number; 
            }
            
            if (debugPrint) System.out.println(Arrays.toString(result));
            br.close();
            
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
