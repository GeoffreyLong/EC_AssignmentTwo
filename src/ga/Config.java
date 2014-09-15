package ga;

import java.awt.Point;
import java.util.List;
import java.util.Map;

import ttp.TTPInstance;
import ttp.TTPSolution;
import ga.Selection.SelectionType;

public class Config{
	public int individualLength = -1;
	public SelectionType selectionType = null;
	public SelectionType parentSelectionType = null;
	
	public int populationSize = 100;
	public int maxNumberOfGenerations = 10000;
	public double crossoverChance = 1;
	public double crossingChance = 1;
	public double mutationChance = 0.02;
	public int tournamentSize = 5;
	public double inverOverProbability = 0.02;
	public boolean generationMix = true;
	
	public double[] mutationTypeChance = new double[5];
	public double[] crossoverTypeChance = new double[4];
	public double[] selectionTypeChance = new double[4];
	
	public TTPInstance ttpInstance = null;
	public double[] tskpW = null;
	private TTPInstance ttp;
	
	public void setMutationTypeChance(double insert, double swap, double invert, double scramble, double inverOver){
		double probability = 0;
		
		mutationTypeChance[0] = (probability+=insert);
		mutationTypeChance[1] = (probability+=swap);
		mutationTypeChance[2] = (probability+=invert);
		mutationTypeChance[3] = (probability+=scramble);
		mutationTypeChance[4] = (probability+=inverOver);
		if (Math.abs(1 - probability) >= 0.00001){
			throw new IllegalArgumentException("Mutation Chances arguments must add up to one");
		}
	}
	
	public double[] getMutationTypeChance(){
		return this.mutationTypeChance;
	}
	
	public void setCrossoverTypeChance(double order, double pmx, double cycle, double edge){
		double probability = 0;
		
		crossoverTypeChance[0] = (probability+=order);
		crossoverTypeChance[1] = (probability+=pmx);
		crossoverTypeChance[2] = (probability+=cycle);
		crossoverTypeChance[3] = (probability+=edge);
		if (Math.abs(1 - probability) >= 0.00001 && probability != 0){
			throw new IllegalArgumentException("Crossover Chances arguments must add up to one");
		}
	}
	
	public double[] getCrossoverTypeChance(){
		return this.crossoverTypeChance;
	}
	
	public void setSelectionTypeChance(double roulette, double tournament, double sus, double elitism){
		double probability = 0;
		
		selectionTypeChance[0] = (probability+=roulette);
		selectionTypeChance[1] = (probability+=tournament);
		selectionTypeChance[2] = (probability+=sus);
		selectionTypeChance[3] = (probability+=elitism);
		if (Math.abs(1 - probability) >= 0.00001){
			throw new IllegalArgumentException("Selection Chances arguments must add up to one");
		}
	}
	
	public double[] getSelectionTypeChance(){
		return this.selectionTypeChance;
	}
	
	private static Config instance = null;
	public static Config getInstance(){
		if(instance == null){
			instance = new Config();
		}
		return instance;
	}
	
	private void Config(){}
	
	public void setNumberOfGenerations(int generations){
		this.maxNumberOfGenerations=  generations;
	}
	
	public int getNumberOfGenerations(){
		return this.maxNumberOfGenerations;
	}
	
	public void setPopulationSize(int size){
		populationSize = size;
	}
	
	public int getPopulationSize(){
		return populationSize;
	}
	
	public void setIndividualLength(int length) {
		individualLength = length;
	}

	public int getIndividualLength(){
		return individualLength;
	}
	void setGenerationMix(boolean mixGenerations){
		this.generationMix = mixGenerations;
	}
	public void setTournamentSize(int tournamentSize){
		this.tournamentSize = tournamentSize;
	}
	
	public void setSelectionType(SelectionType type){
		selectionType = type;
	}
	public SelectionType getSelectionType(){
		return this.selectionType;
	}
	public void setMutationChance(double chance){
		mutationChance = chance;
	}
	public void setCrossoverChance(double chance){
		crossoverChance = chance;
	}
	public void setInverOverProbability(double chance){
		this.inverOverProbability = chance;
	}
	public void setParentSelectionType(SelectionType parentSelectionType){
		this.parentSelectionType = parentSelectionType;
	}
	public SelectionType getParentSelectionType(){
		return this.parentSelectionType;
	}

	
	public void setTTPInstance(TTPInstance ttp){
		this.ttp = ttp;
	}
	public void setTtpInstance(TTPInstance inst) {
		this.ttpInstance = inst;
	}
	public void setTSKPw(double[] d) {
		this.tskpW= d;
	}
	public double calculateFitness(Individual individual){
		/*
		double fitness = 0;
		List<Object> chromosomes = individual.genotype;
		
		String lastChromosome = (String) chromosomes.get(chromosomes.size() - 1);
		Point lastPoint = (Point) possibleAlleles.get(lastChromosome);
		Point currentPoint;
		
		for (Object chromosome : individual.genotype){
			chromosome = (String) chromosome;
			currentPoint = (Point) possibleAlleles.get(chromosome);
			fitness += currentPoint.distance(lastPoint);
			lastPoint = currentPoint;
		}
		*/
		double obj = 0;
		double totalWeight = 0;
		double weightConst = (ttpInstance.maxSpeed - ttpInstance.minSpeed) / ttpInstance.capacityOfKnapsack;
		int cityStart = 0;
		obj += this.ttpInstance.distances(cityStart,Integer.parseInt((String)individual.genotype.get(0))) / ttpInstance.maxSpeed;
		for (int i = 0; i < individual.genotype.size() - 1; i++) {
			int cityA = Integer.parseInt((String) individual.genotype.get(i));
			int cityB = Integer.parseInt((String)individual.genotype.get(i+1));
			totalWeight += tskpW[cityA]; // pick up items at cityA
			obj += ttpInstance.distances(cityA, cityB) / (ttpInstance.maxSpeed - totalWeight*weightConst);
		}
		int cityEnd = Integer.parseInt((String)individual.genotype.get(individual.genotype.size()-1));
		totalWeight += tskpW[cityEnd];
		obj += ttpInstance.distances(cityEnd, cityStart) / (ttpInstance.maxSpeed - totalWeight*weightConst);
		return -ttpInstance.rentingRatio*obj;
	}
	public double calculateFitnessTwo(ttp.newrep.Individual individual){
		TTPSolution solution = new TTPSolution(ttp.getTour(individual), ttp.getPackingPlan(individual));
		ttp.evaluate(solution);
		return solution.ob;
	}
}
