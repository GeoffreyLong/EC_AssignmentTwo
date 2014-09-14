package ga;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;


public class Selection {	
	private SelectionType selectionType;
	private static Random rand = new Random(System.currentTimeMillis());
	
	public enum SelectionType{
		ROULETTE,TOURNAMENT,SUS,ELITISM
	}
	
	public Selection(SelectionType selectionType){
		this.selectionType = selectionType;
	}
	
	public Population select(Population population){
		Config config = Config.getInstance();
		if (selectionType == null){
			return population;
		}
		switch(selectionType){
			case ROULETTE:
				return rouletteWheel(population, config.populationSize);
			case SUS:
				return stochasticUniversalSampling(population, config.populationSize);
			case TOURNAMENT:
				return tournamentSelection(population, config.populationSize, config.tournamentSize);		
			case ELITISM:
				return elitism(population, config.populationSize);
			default:
				return population;
		}
	}
    
    public Population rouletteWheel(Population pop, int outSize){
    	Individual [] subset = new Individual[outSize];
    	double[] maxFitScores = new double[pop.population.size()];
    	
    	double populationFitness = pop.calculateTotalFitness();
    	
    	for (int i = 0; i<pop.population.size(); i++){
    		maxFitScores[i]=Config.getInstance().calculateFitness(pop.population.get(i))/populationFitness;
    	}
    	
    	int outCount=0;
    	while (outCount<outSize){
    		double index = rand.nextDouble()*populationFitness;
    		for (int i = 0; i<pop.population.size(); i++){
    			if(index<=maxFitScores[i]){
    				subset[outCount]=pop.population.get(i);
    				outCount++;
    				break;
    			}
    		}
    	}
    	return new Population(subset);
    }
    
    public Population stochasticUniversalSampling(Population pop, int outSize){
    	Individual [] subset = new Individual[outSize];
    	double[] maxFitScores = new double[pop.population.size()];
    	double populationFitness = pop.calculateTotalFitness();
    	
    	for (int i = 0; i<pop.population.size(); i++){//calculate the max fitness proportion space for each individual
    		if(i==0){
    			maxFitScores[i]=Config.getInstance().calculateFitness(pop.population.get(i))/populationFitness;
    		}else{
    			maxFitScores[i]=maxFitScores[i-1]+Config.getInstance().calculateFitness(pop.population.get(i))/populationFitness;
    		}    		
    	}
    	
    	double index = rand.nextDouble()*(1.0/outSize);//get random start index between 0 and 1/outSize
    	int i = 0;
    	int outCount=0;
    	while (i<pop.population.size()){
    		if (index<=maxFitScores[i]){
    			subset[outCount]=pop.population.get(i);
    			outCount++;
    			index+=1.0/outSize;
    		}else{
    			i++;
    		}
    	}
    	
    	return new Population(subset);
    }
    
    public Population tournamentSelection(Population pop, int outSize, int tourSize){
    	Individual [] subset = new Individual[outSize];
    	int popSize = pop.size();
    	int [] indexes = new int[tourSize];
    	Set<Integer> indexesB = new HashSet<Integer>();    	
    	
    	int outCount = 0;
    	while (outCount<outSize){//until we have the output subset population size
    		int tourCount=0;
    		double bestFitness = 0;
    		int bestIndex = -1;
    		while (tourCount<tourSize){//until we have the specified tour size
    			int index = rand.nextInt(popSize);
    			if (!indexesB.contains(index)){			// If index hasn't yet been selected add it
    				indexesB.add(index);// only need if using probability
    				indexes[tourCount]=index;
    				tourCount++;
    				double fitness=Config.getInstance().calculateFitness(pop.population.get(index));
    				if(fitness>bestFitness){// fitness of this individual is best
    					bestFitness=fitness;
    					bestIndex=index;
    				}// take best one? or take best one with probability prob (as wiki suggests)
    			}
    		}
    		subset[outCount]=pop.population.get(bestIndex);
    		outCount++;
    		indexesB.clear();
    	} 
    	return new Population(subset);
    }
    
    public Population elitism(Population pop, int outSize){//cut percent (rather than number)
    	Individual [] subset = new Individual[outSize];
    	//sort by fitness
    	Comparator<Individual> indComp = new Comparator<Individual>() {
    		@Override
    		public int compare(Individual a, Individual b) {
    			double diff = Config.getInstance().calculateFitness(a) - Config.getInstance().calculateFitness(b);
    			return (diff == 0) ? 0 : ((diff > 0) ? -1 : 1);
    		}
    	};
    	

    	Collections.sort(pop.population,indComp);//java 7
    	//pop.population.sort(indComp);//java 8

    	
    	//cut off
    	for(int i=0; i<outSize; i++){
    		subset[i]=pop.population.get(i);
    	}
    	
    	return new Population(subset);
    }
}
