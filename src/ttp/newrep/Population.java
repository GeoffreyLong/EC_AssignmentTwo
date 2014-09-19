package ttp.newrep;

import java.util.LinkedList;
import java.util.List;

import ttp.TTPInstance;
import ttp.newrep.helpers.*;

public class Population {
	public enum MutationType{ONE, TWO, NONE};
	public enum CrossoverType{ONE, TWO, NONE};
	public enum Fixer{REMOVE_BAD_PROFITS, ONE, TWO, NONE};
	public enum Sorter{ONE, TWO, NONE};
	public enum Selector{ONE, TWO, NONE};
	public enum Packer{GREEDY_PROFIT, GREEDY_PROFIT_SLOW, ONE, TWO, NONE};
	
	public Individual[] population;
	public TTPInstance instance;
	public Population(TTPInstance instance, int[] tour, int popSize){
		for(int i = 0; i < popSize; i++){
			population[i] = instance.createIndividual(tour);
		}
		this.instance = instance;
	}
	
	public void evolve(int maxNumberOfIterations, Sorter sorter, MutationType mutationType, CrossoverType crossoverType, Fixer fixer, Packer packer, Selector selector){
		int iterationNumber = 0;
		while(iterationNumber < maxNumberOfIterations){
			iterationNumber ++;
			List<Individual> newPop = new LinkedList<Individual>();
			for (Individual ind : population){
				newPop.add(ind);
				newPop.add(instance.createIndividual(instance.getTour(ind), instance.getTour(ind)));
			}
			
			switch(sorter){
			
			}
			
			switch(mutationType){
				case ONE:
					newPop = Mutators.mutateSwitchNeighbor(instance, newPop);
					break;
				case TWO:
				case NONE:
			}
			
			switch(crossoverType){
				case ONE:
					
				case TWO:
				case NONE:	
			}
			
			switch(fixer){
				case REMOVE_BAD_PROFITS:
					newPop = Fixers.removeBadProfits(instance, newPop);
					break;
			}
			
			switch(packer){
				case GREEDY_PROFIT:
					newPop = Packers.greedyProfits(instance, newPop);
					break;
				case GREEDY_PROFIT_SLOW:
					newPop = Packers.greedyProfitsSlow(instance, newPop);
					break;
			}
			
			switch(selector){
			
			}
		}
	}
}
