

import java.io.*;

import ttp.Optimisation.Optimisation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import ttp.TTPInstance;
import ttp.TTPSolution;
import ttp.Utils.DeepCopy;
import ttp.Utils.Utils;
import ttp.newrep.Individual;

/**
 *
 * @author wagner
 */
public class Driver {
    
    /* The current sequence of parameters is
     * args[0]  folder with TTP files
     * args[1]  pattern to identify the TTP problems that should be solved
     * args[2]  optimisation approach chosen
     * args[3]  stopping criterion: number of evaluations without improvement
     * args[4]  stopping criterion: time in milliseconds (e.g., 60000 equals 1 minute)
     */
    public static void main(String[] args) {
       
        if (args.length==0) 
        	args = new String[]{"instances", "a280_n279_bounded-strongly-corr_01",
        	//args = new String[]{"instances", "a280_n1395_uncorr-similar-weights_05",
        	//args = new String[]{"instances", "a280_n2790_uncorr_10",
        	//args = new String[]{"instances", "fnl4461_n4460_bounded-strongly-corr_01",
        	//args = new String[]{"instances", "fnl4461_n22300_uncorr-similar-weights_05",
        	//args = new String[]{"instances", "fnl4461_n44600_uncorr_10",
        	//args = new String[]{"instances", "pla33810_n33809_bounded-strongly-corr_01",
        	//args = new String[]{"instances", "pla33810_n169045_uncorr-similar-weights_05",
        	//args = new String[]{"instances", "pla33810_n338090_uncorr_10",
            "11", "5", "600000"};
//        ttp.Optimisation.Optimisation.doAllLinkernTours();
//        runSomeTests();
        doBatch(args);
        //testAllInst();
    }
    
    // note: doBatch can process several files sequentially
    public static void doBatch(String[] args) {;
    	
    	File[] files = ttp.Utils.Utils.getFileList(args);
        
        int algorithm = Integer.parseInt(args[2]);
        int durationWithoutImprovement = Integer.parseInt(args[3]);
        int maxRuntime = Integer.parseInt(args[4]);
        

        
        for (File f:files) {
            // read the TSP instance
            TTPInstance instance = new TTPInstance(f);   
            
            // generate a Linkern tour (or read it if it already exists)
            int[] tour = Optimisation.linkernTour(instance);
            
            long startTime = System.currentTimeMillis();
            String resultTitle="";  

            TTPSolution newSolution=null;
            switch (algorithm){
            	case 1: //Hill-Climber (Provided)
            		newSolution=Optimisation.hillClimber(instance, tour, 2, durationWithoutImprovement, maxRuntime);
            		resultTitle = instance.file.getName() + ".hillClimber." + startTime;
            		break;
            	case 2: //Simple Heuristic (From Research)
            		newSolution=Optimisation.simpleHeuristic(instance, tour, maxRuntime);
            		resultTitle = instance.file.getName() + ".simpleHeuristic." + startTime;
            		break;
            	case 3: //Exercise 2 : Algorithm 1 : Greedy Heuristic Packing Plan Change            		
            		newSolution=Optimisation.ppGreedyRegardTour(instance, tour);            		
            		resultTitle = instance.file.getName() + ".ppGreedyRegardTour." + startTime;
            		break;
            	case 10: //Exercise 2 : Algorithm 1 : Greedy Heuristic Packing Plan Change            		
            		TTPSolution temp=Optimisation.ppGreedyRegardTour(instance, tour);
            		newSolution=Optimisation.bitFlip(instance, temp,(int)(maxRuntime-(System.currentTimeMillis()-startTime+10)),0);            		
            		resultTitle = instance.file.getName() + "ppGreedyRegardTour.bitFlip." + startTime;
            		break;
            	case 4: //Exercise 3 : Algorithm 1 : Greedy Heuristic Packing Plan with Tour Flip Potential 
            		newSolution=Optimisation.flipTourCheck(instance,tour);//check whether should flip and apply original PPlan
            		resultTitle = instance.file.getName() + ".ppGreedyRegardTour_flip." + startTime;
            		break;
            	case 5: //Exercise 3 : Algorithm 2 : Continuous tour building, flipping, and PP assignment
            		newSolution=Optimisation.randomLinkernTours(instance, maxRuntime);
            		resultTitle = instance.file.getName() + ".randomLinkernTours_ppGreedyRegardTour_flip." + startTime;
            		break;
            	case 6: //Exercise 3 : Algorithm 3 : A2 + insertion
            		TTPSolution tmpSolution=Optimisation.randomLinkernTours(instance, maxRuntime-(int)(instance.numberOfItems*.001));
            		newSolution = Optimisation.insertion(instance, tmpSolution.tspTour, tmpSolution.packingPlan, (int)(instance.numberOfItems*.001));
            		resultTitle = instance.file.getName() + ".randomLinkernTours_ppGreedyRegardTour_flip_insert." + startTime;
            		break;
            	case 11://EA VERSION
                    int gen =1;
                    int MAX_GENS=20;
                    int POP_SIZE=100;
                    boolean diffTours=true;//EX 2 = false, EX 3 = true
                    
                    TTPSolution[] population = Optimisation.instantiatePop(instance,tour,POP_SIZE,diffTours,(int)(maxRuntime/5.0));//use linkern
                    
                    TTPSolution[] newPopulation = new TTPSolution[population.length];
                    while (gen<=MAX_GENS && (System.currentTimeMillis() - startTime)<maxRuntime){
                    	Random rand = new Random();
                    	//clone solutions
                    	for(int i=0; i<population.length; i++){
                    		newPopulation[i]=new TTPSolution(population[i].tspTour,population[i].packingPlan);
                    	}
                    	if(!diffTours){
                    		//mutate only PP (EX 2)
                    		for(int i=0; i<population.length; i++){
                    			Optimisation.bitFlip(instance, newPopulation[i], (int)(maxRuntime/MAX_GENS*POP_SIZE), rand.nextDouble());
                    		}
                    	}else{
                    		//mutate entire TTP (EX 3)
                    		//for(int i=0; i<population.length; i++){
                    		//	Optimisation.insertion(instance, newPopulation[i].tspTour, newPopulation[i].packingPlan, (int)(maxRuntime/MAX_GENS*POP_SIZE));
                    		//}
                    	}
                    	
                    	//crossover (EX 4)
                    	if(diffTours){
                    	for(int i=0; i<POP_SIZE/2; i++){
                    		TTPSolution[] tmp = Optimisation.crossoverSubTours(instance, newPopulation[i*2], newPopulation[i*2+1]);
                    		newPopulation[i*2]=tmp[0];
                    		newPopulation[i*2+1]=tmp[1];
                    	}
                    	}
                    	//select pop size best from combined both
                    	//merge sort and cut arrays based on OB
                    	TTPSolution[] merged = new TTPSolution[POP_SIZE*2];
                    	for(int i=0; i<POP_SIZE*2; i++){
                    		if(i<POP_SIZE){
                    			merged[i]=population[i];
                    		}else{
                    			merged[i]=newPopulation[i-POP_SIZE];
                    		}
                    		instance.evaluate(merged[i]);
                    	}
                    	Arrays.sort(merged, new Comparator<TTPSolution>() {
                    	    public int compare(TTPSolution ttp1, TTPSolution ttp2) {
                    	        return Double.compare(ttp2.ob,ttp1.ob);
                    	    }
                    	});
                    	
                    	//update population
                    	System.out.println("GEN: "+gen);
                    	for(int i=0; i<POP_SIZE; i++){
                    		population[i]=merged[i];
                    		System.out.println(population[i].ob);
                    	}
                    	gen++;
                    }
                    //take top of the array above as solution
                    newSolution=population[0];
                    break;
            }
                    
            newSolution.computationTime = System.currentTimeMillis() - startTime;
            newSolution.writeResult(resultTitle);
            newSolution.altPrint();
        }
    }
    
    public static void testAllInst(){
    	CSVReader csvReader;
		try {    	
	    	List<String[]> allArgs = new LinkedList<String[]>();
	    	allArgs.add(new String[]{"instances", "a280_n279_bounded-strongly-corr_01"});
	    	allArgs.add(new String[]{"instances", "a280_n1395_uncorr-similar-weights_05"});
	    	allArgs.add(new String[]{"instances", "a280_n2790_uncorr_10"});
	    	allArgs.add(new String[]{"instances", "fnl4461_n4460_bounded-strongly-corr_01"});
	    	allArgs.add(new String[]{"instances", "fnl4461_n22300_uncorr-similar-weights_05"});
	    	allArgs.add(new String[]{"instances", "fnl4461_n44600_uncorr_10"});
	    	allArgs.add(new String[]{"instances", "pla33810_n33809_bounded-strongly-corr_01"});
	    	allArgs.add(new String[]{"instances", "pla33810_n169045_uncorr-similar-weights_05"});
	    	allArgs.add(new String[]{"instances", "pla33810_n338090_uncorr_10"});
	        
	    	int maxRuntime = 600000;
	    	
	    	for(int i = 0; i < allArgs.size(); i++){
	    		File[] files = ttp.Utils.Utils.getFileList(allArgs.get(i));
	    		for (File f : files){
	    			TTPInstance instance = new TTPInstance(f);
	                // generate a Linkern tour (or read it if it already exists)
	                int[] tour = Optimisation.linkernTour(instance);
	    			
	                
	                //TTPSolution solution = Optimisation.exerciseTwoSolutionOne(instance, tour, instance.createIndividual(tour),1); String csvFilename = "Z_ex2sol1iter1.csv";
	                //TTPSolution solution = Optimisation.exerciseTwoSolutionTwo(instance, tour, 2, maxRuntime,false); String csvFilename = "Z_ex2sol2iter1.csv";
	                //TTPSolution solution = Optimisation.exerciseTwoSolutionTwoAlt(instance, tour, 10, maxRuntime,false); String csvFilename = "Z_ex2sol2Aiter1.csv";
	                //TTPSolution solution = Optimisation.exerciseTwoSolutionTwoAltTwo(instance, tour, 2, maxRuntime, false); String csvFilename = "Z_ex2sol2Biter1.csv";
	                TTPSolution solution = Optimisation.exerciseThreeSolutionTwo(instance, tour, 5, maxRuntime); String csvFilename = "Z_ex3sol2Aiter1.csv";
	                //TTPSolution solution = Optimisation.exerciseFourSolutionOne(instance, tour, 10, maxRuntime); String csvFilename = "Z_ex4sol1iter1.csv";
	    			
	                
	    			List<String[]> data = new ArrayList<String[]>();
	    			csvReader = new CSVReader(new FileReader(csvFilename));
	    			
	    			for (Object line : csvReader.readAll()) {
	    	    	    data.add((String[]) line);
	    	    	}
	                
	    			csvReader.close();
	    			
	                String[] oldData = data.get(i);
	                Double oldOb = Double.valueOf(oldData[1]);
	                long oldCompTime = Long.valueOf(oldData[2]);
	                int iterations = Integer.valueOf(oldData[3]);
	                
	                double ob = (solution.ob + oldOb * iterations) / (iterations + 1);
	                long compTime = (solution.computationTime + oldCompTime * iterations) / (iterations + 1);
	                
	    			String[] newData = new String[]{String.valueOf(i+1), String.valueOf(ob), String.valueOf(compTime), String.valueOf(iterations+1)};
	    			System.out.println("Cur Results // " + (i+1) + " : " + String.format("%.0f", solution.ob) + "\t: " + solution.computationTime);
	    			//System.out.println("Tot Results // " + (i+1) + " : " + String.format("%.0f", ob) + "\t: " + compTime);
	    			data.set(i, newData);
	    			
	    			CSVWriter writer = new CSVWriter(new FileWriter(csvFilename));
	    			writer.writeAll(data);
	    	    	writer.close();
	    		}
	    	}
	    	
		} catch (FileNotFoundException e) {
			String fileName = e.getLocalizedMessage().split(" ")[0];
			System.out.println("FILE NOT FOUND, writing file " + fileName);
			try {
				CSVWriter writer = new CSVWriter(new FileWriter(fileName));
				String[] string = new String[]{"0", "0", "0", "0"};
				List<String[]> data = new ArrayList<String[]>();
				for(int i=0; i<9; i++){
					writer.writeNext(string);
				}
				writer.close();
				testAllInst();
			} catch (IOException e1) {
				System.out.println("Could not write file " + fileName);
				e1.printStackTrace();
			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    public static void runSomeTests() {
        //        TTPInstance instance = new TTPInstance(new File("instances/a280_n279_bounded-strongly-corr_1.ttp"));
        //TTPInstance instance = new TTPInstance(new File("instances/a280_n1395_bounded-strongly-corr_1.ttp"));
        TTPInstance instance = new TTPInstance(new File("instances/a280_n279_uncorr_01.ttp"));
//        TTPInstance instance = new TTPInstance(new File("instances/a280_n2790_bounded-strongly-corr_10.ttp"));
//        TTPInstance instance = new TTPInstance(new File("instances/a280_n837_uncorr_9.ttp"));
//        instance.printInstance(false);
        
        int[] tour = new int[instance.numberOfNodes+1];
//        for (int i=0; i<tour.length; i++) tour[i] = i;
//        tour[instance.numberOfNodes]=0;
////        tour = permutation(tour.length);
        
        ttp.Utils.Utils.startTiming();
        tour = Optimisation.linkernTour(instance);
        ttp.Utils.Utils.stopTimingPrint();
        
        
        int[] packingPlan = new int[instance.numberOfItems];
        TTPSolution solution = new TTPSolution(tour, packingPlan);
        instance.evaluate(solution);
        System.out.print("\nLINKERN tour and no pickup: ");
        solution.printFull();
        
        packingPlan = new int[instance.numberOfItems];
        for (int i=0; i<packingPlan.length; i++) packingPlan[i] = 0;
//        for (int i=0; i<packingPlan.length; i++) packingPlan[i] = Math.random()<0.1?1:0;
        packingPlan[0]=1;
//        packingPlan[11]=1;
//        packingPlan[12]=1;
//        packingPlan[packingPlan.length-1]=1;
//        TTPSolution solution = new TTPSolution(tour, packingPlan);
//        instance.evaluate(solution);
//        solution.print();
        solution = new TTPSolution(tour, packingPlan);
        instance.evaluate(solution);
        System.out.print("\nLINKERN tour and only pickup of the first item: ");
        solution.printFull();
        
        int durationWithoutImprovement = 100;
        
        System.out.println("\nOptimiser: hillclimber (flip 1)");
        Optimisation.hillClimber(instance, tour, 1, durationWithoutImprovement, 600).printFull();
        
        System.out.println("\nOptimiser: hillclimber (flip with prob 1/n)");
        Optimisation.hillClimber(instance, tour, 2, durationWithoutImprovement, 600).printFull();
        
        
    }
    
    
}
