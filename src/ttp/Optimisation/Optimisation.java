package ttp.Optimisation;


import ga.Mutation;
import ga.Population;
import ga.Config;

import java.awt.Point;
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
    
	public static TTPSolution cosolver(TTPInstance instance, int[] tour, int maxRuntime) {
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
		
		while (runtime<maxRuntime){
			packingPlanDash = solveKRP(instance,d,W,individual);
			
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
		}
				
		TTPSolution s = new TTPSolution(tourRet, packingPlanRet);
		long duration = ttp.Utils.Utils.stopTiming();
	    s.computationTime = duration;
    	instance.evaluate(s);
        return s;
        
    }
	

	private static int[] solveKRP(TTPInstance instance, double[] d, double[] W, Individual individual){
		
		int[] packingPlanRet = new int[instance.numberOfItems];
		double P = Double.NEGATIVE_INFINITY;
		double profit = 0;
		double t = 0;
		int[][] items= instance.items;		
		int itemsPerCity = instance.numberOfItems / individual.tour.length;
		
		//calc profit
		System.out.println("VALUES MUST BE EQUAL: "+individual.tour[0].items.size()+" : "+itemsPerCity);
		for (int i = 0; i < individual.tour.length; i++){
			for(int j = 0; j < individual.tour[i].items.size(); j++){
				//if(packingPlan[(i*itemsPerCity + j)]==1)
				if(individual.tour[i].items.get(j).isSelected){
					profit += individual.tour[i].items.get(j).profit;
				}
			}
		}
		
		//calc renting rate
		for (int i = 0; i < instance.numberOfNodes; i++){
			t += d[i]/(instance.maxSpeed - W[i]*((instance.maxSpeed-instance.minSpeed)/instance.capacityOfKnapsack));
		}
		
		P = profit - instance.rentingRatio*t;
		
		return packingPlanRet;
	}
	private static int[] solveTSKP(double[] W, ttp.newrep.Individual ind) {
		Config config = Config.getInstance();
		config.setTSKPw(W);
		int populationSize = 100;
		config.setInverOverProbability(0.02);
		// set inverOver probability and fitness function
		Mutation mutate = new Mutation(null);
		Population population = new Population(populationSize-1, ind.tour.length);
		
		ga.Individual currentSol = new ga.Individual();
		currentSol.genotype = new ArrayList<Object>(ind.tour.length);
		for (int i = 0; i < ind.tour.length; i++) {
			currentSol.genotype.set(i, ind.tour[i].cityId);
		}
		population.population.add(currentSol);
		
		
		int numberOfGenerations = 0;
		int maxGeneration = 1000;
		double bestSolution = Double.MAX_VALUE;
		ga.Individual bestSolInd=null;
		
		System.out.println("--------------------------------------------------------------------------");
		System.out.println("GEN #     ITER BEST (  POP BEST,    POP AVG,  POP WORST), TIME SINCE ITER START ( OVERALL TIME AVG, OVERALL TIME SUM)");
		while (numberOfGenerations < maxGeneration){
			Population offspring = population.clone();
			
			offspring = mutate.inverOver(offspring);
			
			numberOfGenerations++;
			
			/// calc data store best worst and avg
			//bestF,avgF,worstF,bestInd,worstInd
			Double[] data = population.getStats();

			if (data[0] < bestSolution){
				bestSolution = data[0];
				bestSolInd=population.population.get(data[3].intValue());
			}
			
			System.out.println("G: "+String.format("%5d",numberOfGenerations)+" "+String.format("%10.3f",bestSolution) + " ("+String.format("%10.2f",data[0])+", "+String.format("%10.2f",data[1])+", "+String.format("%10.2f",data[2])+"), \n");
			
		}
		
		int[] sol = new int[bestSolInd.genotype.size()+1];
		sol[0] = 1;
		for (int i = 1; i < sol.length; i++) {
			sol[i] = (int)bestSolInd.genotype.get(i);
		}
		
		return sol;
}
	
    public static TTPSolution simpleHeuristic(TTPInstance instance, int[] tour, int maxRuntime) {
    	ttp.Utils.Utils.startTiming();
    	double[] D = new double[instance.numberOfNodes];
    	double dSum = 0;
    	D[instance.numberOfNodes-1] = 0; 
    	for (int i = instance.numberOfNodes-2; i >= 0; i--) {
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

    			packingPlan[ppIndex] = 1;
    			
    			Wc += instance.items[itemIdx[i]][2];
    			//System.out.println("i: "+arrIndex+" ppI: "+ppIndex + " CI " + (cityIndex) + " id: " + (instance.items[itemIdx[i]][0])+" IN: "+itemNumber+" WGC: "+ instance.items[itemIdx[i]][2] +" WGT: "+Wc);
    			//System.out.printf("ItemID: %d\n", itemIdx[i]);
    		}
    		if (Wc == instance.capacityOfKnapsack) {
    			break;
    		}
    	}
    	//System.out.println("Our wc: " +Wc);
    	
    	
   
    	//System.out.println(Arrays.toString(packingPlan));
    	TTPSolution s = new TTPSolution(tour, packingPlan);
        long duration = ttp.Utils.Utils.stopTiming();
        s.computationTime = duration;
    	instance.evaluate(s);	
        return s;
    }
    
    public static TTPSolution getCuts(TTPInstance instance, int[] tour){
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
		    		}*/
				}
    		}
    	}
    	int[] packingPlan = new int[instance.numberOfItems];
    	int weight = 0;
    	for (int k = 0; k < items.size(); k++){ 		
			//System.out.println(items.get(k)[1] + " " + items.get(k)[0] + " " + items.get(k)[2]);
			//System.out.println(items.get(k)[0] + " " + tour[k]);
			if (items.get(k)[1] > weight){
				weight += items.get(k)[2];
				packingPlan[k] = 1;
			}
			
			TTPSolution solution = new TTPSolution(tour, packingPlan);
	    	instance.evaluate(solution);
	    	System.out.println(solution.ob + " " + solution.wendUsed);
		}
    	for (int k = 0; k < tour.length; k++){
    		
    	}
    	
    	return null;
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
    public static TTPSolution exerciseTwoSolutionTwo(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
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
		
    	for (int i = tour.length-1; i >= 0; i--){
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
    	
    	int[] packingPlan = new int[instance.numberOfItems];
    	int[] packingPlanClone = packingPlan.clone();
    	TTPSolution newSolution = new TTPSolution(tour, packingPlan);
        instance.evaluate(newSolution);
        double lastSolutionOb = -Double.MAX_VALUE;
        int didNotImprove = 0;
		
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
    		
    		lastSolutionOb = newSolution.ob;
    		System.out.println(newSolution.ob);
    		
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
    
    /**
     * Swaps random two random nodes and mutates the packing plans on these nodes stochastically
     * If the fitness of the new tour and plan is better than the last, this tour/plan will be used in future mutations
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseThreeSolutionOne(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
    	int[] packingPlan = new int[instance.numberOfItems];
    	TTPSolution newSolution = new TTPSolution(tour, packingPlan);
        instance.evaluate(newSolution);
    	
        double lastSolutionOb = -Double.MAX_VALUE;
        int didNotImprove = 0;
		
        while(didNotImprove <= durationWithoutImprovement){
    		if (lastSolutionOb == newSolution.ob){
    		}
    		else{
    			didNotImprove = 0;
    		}
    		
    		lastSolutionOb = newSolution.ob;
    		System.out.println(newSolution.ob);
    		
    		int randNodeIndexOne = (int) (Math.random() * tour.length);
			int randNodeIndexTwo = (int) (Math.random() * tour.length);
			if (randNodeIndexOne != 0 && randNodeIndexOne != tour.length-1 && randNodeIndexTwo != 0 && randNodeIndexTwo != tour.length-1){
				int itemsPerCity = instance.numberOfItems / (tour.length - 2);
				int[] packingListIndices = new int[2*itemsPerCity];
				
	    		for (int i = 0; i < 5; i++){
	    			int[] tourClone = tour.clone();
	    			int[] packingPlanClone = packingPlan.clone();
	    			double randSwap = Math.random();
	    			if (randSwap <= 0.1){
	    				tourClone[randNodeIndexOne] = tour[randNodeIndexTwo];
	    				tourClone[randNodeIndexTwo] = tour[randNodeIndexOne];
	    			}
	    			if (tourClone[randNodeIndexOne] != 0 && tourClone[randNodeIndexTwo] != 0){
		    			for (int j = 0; j < itemsPerCity; j++){
		    				double planMutate = Math.random();
		    				int itemIndex = 0;
		    				if (packingListIndices.length / (j+1) == 2){
		    					itemIndex = (tour.length-2) * j/2 + tourClone[randNodeIndexOne]-1;
		    				}
		    				else{
		    					itemIndex = (tour.length-2) * j + tourClone[randNodeIndexTwo]-1;
		    				}
		    				if (planMutate <= 0.1){
		    					if (packingPlanClone[itemIndex] == 1) packingPlanClone[itemIndex] = 0;
		    					else packingPlanClone[itemIndex] = 1;
		    				}
		    			}
	    			}
	    			TTPSolution tempSol = new TTPSolution(tourClone, packingPlanClone);
	    	        instance.evaluate(tempSol);
	    	        
	    			if (tempSol.ob > lastSolutionOb){
	    				lastSolutionOb = tempSol.ob;
	    				packingPlan = packingPlanClone.clone();
	    				tour = tourClone.clone();
	    			}
	    		}
			}
	    	newSolution = new TTPSolution(tour, packingPlan);
    		instance.evaluate(newSolution);
        }
        
    	return newSolution;
    }
    
    /**
     * Same as the last algorithm, except it uses only adjacent nodes
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseThreeSolutionTwo(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
    	int[] packingPlan = new int[instance.numberOfItems];
    	TTPSolution newSolution = new TTPSolution(tour, packingPlan);
        instance.evaluate(newSolution);
    	
        double lastSolutionOb = -Double.MAX_VALUE;
        int didNotImprove = 0;
		
        while(didNotImprove <= durationWithoutImprovement){
    		if (lastSolutionOb == newSolution.ob){
    		}
    		else{
    			didNotImprove = 0;
    		}
    		
    		lastSolutionOb = newSolution.ob;
    		System.out.println(newSolution.ob);
    		
    		int randNodeIndex = (int) (Math.random() * tour.length);
			if (randNodeIndex > 1 && randNodeIndex < tour.length-1){
				int itemsPerCity = instance.numberOfItems / (tour.length - 2);
				
	    		for (int i = 0; i < 5; i++){
	    			int[] tourClone = tour.clone();
	    			int[] packingPlanClone = packingPlan.clone();
	    			double randSwap = Math.random();
	    			if (randSwap <= 0){
	    				tourClone[randNodeIndex] = tour[randNodeIndex -1];
	    				tourClone[randNodeIndex - 1] = tour[randNodeIndex];
	    			}
    				for (int j = 0; j < itemsPerCity; j++){
	    				double planMutate = Math.random();
	    				int itemIndex = (int)(Math.random() * itemsPerCity) + tourClone[randNodeIndex] - 1;
	    				itemIndex -= (int) (Math.random() + 0.5);
	    				if (planMutate <= 0.5){
	    					if (packingPlanClone[itemIndex] == 1) packingPlanClone[itemIndex] = 0;
	    					else packingPlanClone[itemIndex] = 1;
	    				}
	    			}
 
	    			TTPSolution tempSol = new TTPSolution(tourClone, packingPlanClone);
	    	        instance.evaluate(tempSol);
	    	        
	    			if (tempSol.ob > lastSolutionOb){
	    				lastSolutionOb = tempSol.ob;
	    				packingPlan = packingPlanClone.clone();
	    				tour = tourClone.clone();
	    			}
	    		}
			}
	    	newSolution = new TTPSolution(tour, packingPlan);
    		instance.evaluate(newSolution);
        }
        
    	return newSolution;
    } 
    
    /**
     * Same as last algorithm, but uses the object oriented approach
     * This one does not work out too well
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseThreeSolutionTwoNew(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
        Individual individual = instance.createIndividual(tour);
        int didNotImprove = 0;
        TTPSolution solution = null;
        int[] packingPlan = new int[instance.numberOfItems];
		
        while(didNotImprove <= durationWithoutImprovement){
    		int randNodeIndex = (int) (Math.random() * individual.tour.length);
    		
    		solution = new TTPSolution(tour, packingPlan);
    		instance.evaluate(solution);
    		
			if (randNodeIndex > 1 && randNodeIndex < individual.tour.length){
	    		for (int i = 0; i < 5; i++){
	    			int[] tourNew = tour.clone();
	    			int[] packingPlanNew = packingPlan.clone();
	    			
	    			Individual tempInd = instance.createIndividual(tourNew, packingPlanNew);
	    			double randSwap = Math.random();
	    			if (randSwap <= 0.1){
	    				City cityTemp = tempInd.tour[randNodeIndex - 1];
	    				tempInd.tour[randNodeIndex - 1] = individual.tour[randNodeIndex];
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
	    			if (newSolution.ob > solution.ob && newSolution.wend >= 0){
	    	    		//System.out.println(solution.ob);
	    	    		tour = instance.getTour(tempInd).clone();
	    	    		packingPlan = instance.getPackingPlan(tempInd).clone();
	    	    		didNotImprove = 0;
	    			}
	    			else{
	    				didNotImprove ++;
	    			}
	    		}
			}
        }
        
    	return solution;
    } 
    
    /**
     * Might actually be the same as Two New, but it uses a resetter for the packing plan... 
     * It actually works very well, so definitely keep this one
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseThreeSolutionTwoAlt(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
        Individual individual = instance.createIndividual(tour);
        int didNotImprove = 0;
        TTPSolution solution = null;
        int[] packingPlan = new int[instance.numberOfItems];
        double oldFitness = 0;
        int tempCount = 0;
        double bestFitness = 0;
		
        while(didNotImprove <= durationWithoutImprovement){
    		int randNodeIndex = (int) (Math.random() * individual.tour.length);
    		
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
        		System.out.println(bestFitness);
    			bestFitness = oldFitness;
    		}
    		
			if (randNodeIndex > 1 && randNodeIndex < individual.tour.length){
	    		for (int i = 0; i < 5; i++){
	    			int[] tourNew = tour.clone();
	    			int[] packingPlanNew = packingPlan.clone();
	    			
	    			Individual tempInd = instance.createIndividual(tourNew, packingPlanNew);
	    			double randSwap = Math.random();
	    			if (randSwap <= 0.1){
	    				City cityTemp = tempInd.tour[randNodeIndex - 1];
	    				tempInd.tour[randNodeIndex - 1] = individual.tour[randNodeIndex];
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
	    	        if (tempCount >= 1000){
	    	        	packingPlan = new int[packingPlan.length];
	    	        	tempCount = 0;
	    	        }
	    			if (newSolution.ob > solution.ob && newSolution.wend >= 0){
	    	    		tour = instance.getTour(tempInd).clone();
	    	    		packingPlan = instance.getPackingPlan(tempInd).clone();
	    			}
	    		}
			}
        }
        
    	return solution;
    } 
    
    
    /**
     * Unfinished Algorithm.  Will use crossover to mutate the tours
     * 	Then it will use hill climber to evaluate the new tour
     * 	if the result of this hill climber is greater than the previous, this tour is the new one
     * 	else re-iterate
     * 
     * This version will have 10 individuals per iteration
     * 
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseThreeSolutionThree(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
        Individual individual = instance.createIndividual(tour);
        int didNotImprove = 0;
        TTPSolution solution = null;
        int[] packingPlan = new int[instance.numberOfItems];
        double oldFitness = hillClimber(instance, tour, 2,1000, 60000).ob;
		
        while(didNotImprove <= durationWithoutImprovement){
    		
    		solution = new TTPSolution(tour, packingPlan);
    		instance.evaluate(solution);
    		
    		Individual[] individuals = new Individual[10];
    		for (int i = 0; i < 10; i++){
				int[] tourNew = tour.clone();
				int[] packingPlanNew = packingPlan.clone();
				Individual tempInd = instance.createIndividual(tourNew, packingPlanNew);
				
				
				int randNodeIndex = (int) (Math.random() * individual.tour.length);
				int randNodeIndexTwo = (int) (Math.random() * individual.tour.length);
	
				int j=0;
				while(true) {
					int indexA = (randNodeIndex+1+j) % individual.tour.length;
					int indexB = (randNodeIndexTwo-j + individual.tour.length) % individual.tour.length; 
	
					City temp = tempInd.tour[indexA];//store temp
					tempInd.tour[indexA] = tempInd.tour[indexB];
					tempInd.tour[indexB] = temp;
					if (Math.abs(indexA-indexB) <= 1 || Math.abs(indexA-indexB) >= individual.tour.length-1) break;
					
					j++;
				}
				
				individuals[i] = tempInd;
    		}
    		int bestIndex = -1;
			for (int i = 0; i<individuals.length; i++){
				Individual tempInd = individuals[i];
				TTPSolution hillSol = hillClimber(instance, instance.getTour(tempInd), 2,10000, 60000);
				System.out.println(hillSol.ob + " " + oldFitness);
				if (hillSol.ob > oldFitness){
					oldFitness = hillSol.ob;
					bestIndex = i;
				}
				i++;
			}
			System.out.println();
			if (bestIndex != -1) tour = instance.getTour(individuals[bestIndex]).clone();
        }
        
    	return solution;
    } 
    
    /**
     * Basically the same as ex 3 sol 3
     * @param instance
     * @param tour
     * @param durationWithoutImprovement
     * @param maxRuntime
     * @return
     */
    public static TTPSolution exerciseThreeSolutionThreeAlt(TTPInstance instance, int[] tour, int durationWithoutImprovement, int maxRuntime){
        Individual individual = instance.createIndividual(tour);
        int didNotImprove = 0;
        TTPSolution solution = null;
        int[] packingPlan = new int[instance.numberOfItems];
        double oldFitness = 0;
        int tempCount = 0;
		
        while(didNotImprove <= durationWithoutImprovement){
    		
    		solution = new TTPSolution(tour, packingPlan);
    		instance.evaluate(solution);
    		
			int[] tourNew = tour.clone();
			int[] packingPlanNew = packingPlan.clone();
			Individual tempInd = instance.createIndividual(tourNew, packingPlanNew);
			
			
			
			for (int i = 0; i<2; i++){
				int randNodeIndex = (int) (Math.random() * individual.tour.length);
				int randNodeIndexTwo = (int) (Math.random() * individual.tour.length);
				if (randNodeIndex != randNodeIndexTwo){
					City cityTemp = tempInd.tour[randNodeIndex];
					tempInd.tour[randNodeIndex] = individual.tour[randNodeIndexTwo];
					tempInd.tour[randNodeIndexTwo] = cityTemp;
				}
			}
			
			TTPSolution hillSol = hillClimber(instance, instance.getTour(tempInd), 2,1000, 10000);
			System.out.println(hillSol.ob);
			if (hillSol.ob > oldFitness){
				oldFitness = hillSol.ob;
				tour = instance.getTour(tempInd).clone();
			}
        }
        
    	return solution;
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
