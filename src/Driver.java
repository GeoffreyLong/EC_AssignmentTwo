

import java.io.*;

import ttp.Optimisation.Optimisation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import ttp.TTPInstance;
import ttp.TTPSolution;
import ttp.Utils.DeepCopy;
import ttp.Utils.Utils;

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
            "2", "10000", "60000"};
//        ttp.Optimisation.Optimisation.doAllLinkernTours();
//        runSomeTests();
        doBatch(args);
        
        //testAllInst();
    }
    
    // note: doBatch can process several files sequentially
    public static void doBatch(String[] args) {
//        String[] args = new String[]{"instances/","a2"};                      // first argument: folder with TTP and TSP files, second argument: partial filename of the instances to be solved   
//        System.out.println("parameters: "+Arrays.toString(args));
    	
    	File[] files = ttp.Utils.Utils.getFileList(args);
        
        int algorithm = Integer.parseInt(args[2]);
        int durationWithoutImprovement = Integer.parseInt(args[3]);
        int maxRuntime = Integer.parseInt(args[4]);
        
//        System.out.println("files.length="+files.length+" algorithm="+algorithm+" durationWithoutImprovement="+durationWithoutImprovement);
//        System.out.println("wend wendUsed fp ftraw ft ob computationTime");
        
        for (File f:files) {
            // read the TSP instance
            TTPInstance instance = new TTPInstance(f);
            
            long startTime = System.currentTimeMillis();
            String resultTitle="";
            
            // generate a Linkern tour (or read it if it already exists)
            int[] tour = Optimisation.linkernTour(instance);

            
            System.out.println(f.getName()+": ");
            
            // do the optimisation
            
            //System.out.println("HILL CLIMBER SOLUTION --------------------------------------");
            //TTPSolution solution = Optimisation.hillClimber(instance, tour, algorithm,durationWithoutImprovement, maxRuntime);
           
           //solution.println();
         	//solution.altPrint();
           //solution.printFull();

            //solution2.altPrint();

            // print to file
            //resultTitle = instance.file.getName() + ".SimpleHeuristic." + startTime;
            //solution2.writeResult(resultTitle);
            
            /*
            System.out.println("E2-A1 (Hayden's) -------------------------------------");
            //TTPSolution solution3 = Optimisation.exerciseTwoSolutionOne(instance, tour, instance.createIndividual(tour),1);
            //TTPSolution solution3 = Optimisation.cosolver(instance, tour, maxRuntime);
            //TTPSolution solution3 = Optimisation.exerciseThreeSolutionH(instance, tour, maxRuntime,1);
            TTPSolution solution3 = Optimisation.backFourth(instance,tour, 60000,1);
            resultTitle = instance.file.getName() + ".exerciseThreeSolutionH." + startTime;
            solution3.writeResult(resultTitle);
            solution3.altPrint();
            //solution2.printFull();
            */

            
            
            //solution.altPrint();


            //TTPSolution solution3 = Optimisation.exerciseTwoSolutionTwo(instance, tour, 5, maxRuntime);
            //solution3.printFull();
            //solution3.altPrint();
            //TTPSolution solution3 = Optimisation.exerciseTwoSolutionTwo(instance, tour, 10, maxRuntime);
            //solution3.printFull();
            //solution3.altPrint();
            //resultTitle = instance.file.getName() + ".exerciseTwoSolutionTwo." + startTime;
            //solution3.writeResult(resultTitle);
            
            
            
            // General form of result output
            // test	: obAVG	: timeAVG :	# of times run
            
            
            //TTPSolution solution = Optimisation.exerciseTwoSolutionOne(instance, tour, instance.createIndividual(tour),1);
            //TTPSolution solution = Optimisation.exerciseTwoSolutionTwo(instance, tour, 2, maxRuntime, true);
            //TTPSolution solution = Optimisation.exerciseTwoSolutionTwoAlt(instance, tour, 10, maxRuntime,true);
            //TTPSolution solution = Optimisation.exerciseThreeSolutionOne(instance, tour, 10, maxRuntime);
            //TTPSolution solution = Optimisation.exerciseThreeSolutionTwo(instance, tour, 10, maxRuntime);
            //TTPSolution solution = Optimisation.exerciseThreeSolutionTwoNew(instance, tour, 10, maxRuntime);
            //TTPSolution solution = Optimisation.exerciseThreeSolutionTwoAlt(instance, tour, 10, maxRuntime);
            //TTPSolution solution = Optimisation.exerciseThreeSolutionThree(instance, tour, 10, maxRuntime);
            //TTPSolution solution = Optimisation.exerciseThreeSolutionH(instance, tour, 10, maxRuntime);
            //TTPSolution solution = Optimisation.exerciseFourSolutionOne(instance, tour, 10, maxRuntime);
            
            //solution.altPrint();
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
	                TTPSolution solution = Optimisation.exerciseTwoSolutionTwoAlt(instance, tour, 10, maxRuntime,false); String csvFilename = "Z_ex2sol2Aiter1.csv";
	                //TTPSolution solution = Optimisation.exerciseThreeSolutionOne(instance, tour, 10, maxRuntime); String csvFilename = "Z_ex3sol1iter1.csv";
	                //TTPSolution solution = Optimisation.exerciseThreeSolutionTwo(instance, tour, 10, maxRuntime); String csvFilename = "Z_ex3sol2iter1.csv";
	                //TTPSolution solution = Optimisation.exerciseThreeSolutionTwoNew(instance, tour, 10, maxRuntime); String csvFilename = "Z_ex3sol2Niter1.csv";
	                //TTPSolution solution = Optimisation.exerciseThreeSolutionTwoAlt(instance, tour, 10, maxRuntime); String csvFilename = "Z_ex3sol2Aiter1.csv";
	                //TTPSolution solution = Optimisation.exerciseThreeSolutionThree(instance, tour, 10, maxRuntime); String csvFilename = "Z_ex3sol3iter1.csv";
	                //TTPSolution solution = Optimisation.exerciseThreeSolutionH(instance, tour, 10, maxRuntime); String csvFilename = "Z_ex3solHiter1.csv";
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
