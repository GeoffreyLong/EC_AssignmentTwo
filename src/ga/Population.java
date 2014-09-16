package ga;

import java.util.ArrayList;
import java.util.List;

public class Population {
	public List<Individual> population;
	
	public Population(){
		population = new ArrayList<Individual>();
	}
	
	public Population(List<Individual> individuals){
		population = individuals;
	}
	
	public Population(Individual[] individualArray){
		population = new ArrayList<Individual>();
		for (Individual i : individualArray){
			population.add(i);
		}
	}
	
	public Population(int popSize){
		population = new ArrayList<Individual>();
		for (int i = 0; i < popSize; i++){
			population.add(new Individual());
		}
	}
	
	public Population(int popSize, int indLength){
		population = new ArrayList<Individual>();
		for (int i = 0; i < popSize; i++){
			population.add(new Individual(indLength,indLength));
		}
	}
	
	public int size(){
		return population.size();
	}
	
	public String toString(){
		String temp = "";
		for (int i = 0; i < population.size(); i++){
			temp += population.get(i).toString();
			temp += System.getProperty("line.separator");
		}
		return temp;
	}
	
	public double calculateMeanFitness(){
		return calculateTotalFitness() / population.size();
	}
	
	public double calculateTotalFitness(){
		double totalFitness = 0;
		for (Individual i : population){
			double indFitness = Config.getInstance().calculateFitness(i);
			totalFitness += indFitness;
		}
		return totalFitness;
	}
	
	public Individual getBest(){
		
		double bestFitness = Double.MIN_VALUE;
		Individual bestI=null;
		for (Individual i : population){
			double indFitness = Config.getInstance().calculateFitness(i);
			
			if(indFitness>bestFitness){
				bestI = i;
			}
		}
		return bestI;
	}
	
	public Individual getWorst(){
		
		double worstFitness = Double.MAX_VALUE;
		Individual worstI=null;
		for (Individual i : population){
			double indFitness = Config.getInstance().calculateFitness(i);
			
			if(indFitness>worstFitness){
				worstI = i;
			}
		}
		return worstI;
	}
	
	public Double[] getStats(){//might be ugly but only loop through population once to get all the stats
		double bestFitness = Double.NEGATIVE_INFINITY;
		double worstFitness = Double.POSITIVE_INFINITY;
		double averageFitness = 0;
		int bestIndex = -1;
		int worstIndex = -1;
		
		for (int i = 0; i<population.size(); i++){
			double indFitness = Config.getInstance().calculateFitness(population.get(i));
			
			if(indFitness<worstFitness){
				
				worstIndex = i;
				worstFitness = indFitness;
			}
			
			if(indFitness>bestFitness){
				bestIndex = i;
				bestFitness = indFitness;
			}
			
			averageFitness+=indFitness;
		}
		
		averageFitness=averageFitness/population.size();
		
		Double[] stats = {1/bestFitness,1/averageFitness,1/worstFitness,(double) bestIndex, (double) worstIndex};
		return stats;
	}
	
	@Override
	public Population clone() {
		List<Individual> newPopList = new ArrayList<Individual>();
		for (Individual i : this.population) {
			newPopList.add(i.clone());
		}
		return new Population(newPopList);
	}
}
