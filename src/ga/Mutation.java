package ga;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

public class Mutation {
	private static Random rand = new Random(System.currentTimeMillis());
	private double[] mutationTypeChance;
	
	public enum MutationType{
		INSERT,SWAP,INVERSION,SCRAMBLE,INVEROVER, INSERT_SPECIAL
	}
	
	public Mutation(double[] mutationTypeChance){
		this.mutationTypeChance = mutationTypeChance;
	}
	
	public Population mutate(Population population){
		Config config = Config.getInstance();
		MutationType mutationType = null;
		double random = rand.nextDouble();
		
		for (int i = 0; i < MutationType.values().length; i++){
			if (random < mutationTypeChance[i]){
				mutationType = MutationType.values()[i];
				break;
			}
		}

		switch(mutationType){
			case INSERT:
				for (int i = 0; i<population.size();i++){
					if (rand.nextDouble() < config.mutationChance){
						population.population.set(i, insert(population.population.get(i)));
					}
				}
				break;
			case SWAP:
				for (int i = 0; i<population.size();i++){
					if (rand.nextDouble() < config.mutationChance){
						population.population.set(i, swap(population.population.get(i)));
					}
				}
				break;		
			case INVERSION:
				for (int i = 0; i<population.size();i++){					
					if (rand.nextDouble() < config.mutationChance){
						population.population.set(i, inversion(population.population.get(i)));
					}
				}
				break;
			case SCRAMBLE:
				for (int i = 0; i<population.size();i++){
					if (rand.nextDouble() < config.mutationChance){
						population.population.set(i, scramble(population.population.get(i)));
					}
				}
				break;
			case INVEROVER:
				population = inverOver(population);
				break;
			case INSERT_SPECIAL:
				for (int i = 0; i<population.size();i++){
					if (rand.nextDouble() < config.mutationChance){
						population.population.set(i, insertSpecial(population.population.get(i)));
					}
				}
				break;
		}
		return population;
	}
	
	public Individual insert(Individual i){
		int numChromosomes=i.genotype.size();
		
		int indexA = rand.nextInt(numChromosomes);
		int indexB = rand.nextInt(numChromosomes);
		
		if (indexA > indexB) {//make sure indexes in ascending order
			int tmp = indexA;
			indexA = indexB;
			indexB = tmp;
		}
		//System.out.println(indexA+", "+indexB);//testing
		
		for (int j = indexA+1; j < indexB; j++) {
			Object temp = i.genotype.get(j);
			i.genotype.set(j, i.genotype.get(indexB));
			i.genotype.set(indexB,temp);			
		}

		return i;
	}
	
	public Individual insertSpecial(Individual i){
		Config config = Config.getInstance();
		double totalWeight = 0;
		double[] weightRatio = new double[config.tskpW.length];
		for (int j = 0; j < config.tskpW.length; j++) {
			totalWeight += config.tskpW[j];
		}
		for (int j = 0; j < config.tskpW.length; j++) {
			weightRatio[j] = config.tskpW[j] / totalWeight;
		}
		
		
		int numChromosomes=i.genotype.size();
		
		int indexA = rand.nextInt(numChromosomes);
		int indexB = (int)Math.floor(weightRatio[indexA]*rand.nextDouble()*numChromosomes);
		
		if (indexA < indexB) {//make sure indexes in ascending order
			for (int j = indexA+1; j < indexB; j++) {
				Object temp = i.genotype.get(j);
				i.genotype.set(j, i.genotype.get(indexB));
				i.genotype.set(indexB,temp);			
			}
		} else {
			for (int j = indexB-1; j > indexA; j--) {
				Object temp = i.genotype.get(j);
				i.genotype.set(j, i.genotype.get(indexA));
				i.genotype.set(indexA,temp);			
			}
		}
		//System.out.println(indexA+", "+indexB);//testing
		
		

		return i;
	}
	
	
	
	public Individual swap(Individual i){
		int numChromosomes=i.genotype.size();
		
		int indexA = rand.nextInt(numChromosomes);
		int indexB = rand.nextInt(numChromosomes);
		//System.out.println(indexA+", "+indexB);//testing
		
		Object temp = i.genotype.get(indexA);//store temp
		i.genotype.set(indexA,i.genotype.get(indexB));
		i.genotype.set(indexB,temp);
		
		return i;	
	}
	
	public Individual inversion(Individual i){
		int numChromosomes=i.genotype.size();
		
		int indexA = rand.nextInt(numChromosomes);
		int indexB = rand.nextInt(numChromosomes);
		if (indexA > indexB) {//make sure indexes in ascending order
			int tmp = indexA;
			indexA = indexB;
			indexB = tmp;
		}
		int swaps = (int) (Math.floor(indexB-indexA+1)/2);//how many swap operations
		//System.out.println(indexA+", "+indexB+", "+swaps);//testing
		
		for (int j = 0; j < swaps; j++) {
			Object temp = i.genotype.get(indexA+j);//store temp
			i.genotype.set(indexA+j,i.genotype.get(indexB-j));
			i.genotype.set(indexB-j,temp);
		}	
		
		return i;	
	}
	
	public Individual scramble(Individual i){
		int numChromosomes=i.genotype.size();
		
		int numberOfScrambles = rand.nextInt(numChromosomes);				// Number of chromosomes to scramble
		int [] indexes = new int[numberOfScrambles];						// Random array of indexes to scramble
		int [] sortedIndexes = new int[numberOfScrambles];					// The same indexes sorted
		Set<Integer> indexesB = new HashSet<Integer>();	// HashSet of the indexes to ensure all are unique
		
		//System.out.println(numberOfScrambles);//testing
		
		int count=0;
		while (count < numberOfScrambles){
			int index = rand.nextInt(numChromosomes);
			if (!indexesB.contains(index)){			// If index hasn't yet been selected add it
				indexesB.add(index);
				indexes[count]=index;
				sortedIndexes[count]=index;
				count++;
			}
		}
		
		Arrays.sort(sortedIndexes);
		
		// Store temp chromosomes
		Map<Integer, Object> temp = new HashMap<Integer, Object>();
		for (int j = 0; j < numberOfScrambles; j++){
			temp.put(sortedIndexes[j],i.genotype.get(sortedIndexes[j]));
			//System.out.print(sortedIndexes[j]+",");//testing
		}
		//System.out.println();//testing
		// Fill in new individual
		for (int j = 0; j < numberOfScrambles; j++){
			i.genotype.set(indexes[j], temp.get(sortedIndexes[j]));
		}
		
		return i;	
	}
	
	public Population inverOver(Population p){

		for (int individualIndex = 0; individualIndex < p.size(); individualIndex++){
			Individual originalIndividual = p.population.get(individualIndex);
			Individual clonedIndividual = originalIndividual.clone();
			
			int genotypeSize = clonedIndividual.genotype.size();
			int maxIndex = genotypeSize - 1;
			
			int firstCityIndex = rand.nextInt(genotypeSize);
			Object firstCity = clonedIndividual.genotype.get(firstCityIndex);
			
			while(true){
				Object secondCity = firstCity;
				
				if(rand.nextDouble() <= Config.getInstance().inverOverProbability){
					while(secondCity.equals(firstCity)){
						secondCity = clonedIndividual.genotype.get(rand.nextInt(genotypeSize));
					}
				}
				else{
					Individual randomIndividual = p.population.get(rand.nextInt(p.size()));
					int index = getIndexOfElement(firstCity, randomIndividual.genotype);
					secondCity = randomIndividual.genotype.get((index + 1) % genotypeSize);
				}
				
				int secondCityIndex = getIndexOfElement(secondCity, clonedIndividual.genotype);
				int indexDifference = Math.abs(firstCityIndex - secondCityIndex);
				if (indexDifference == 1 || indexDifference == maxIndex){
					break;
				}
				
				int j=0;
				while(true) {
					int indexA = (firstCityIndex+1+j) % genotypeSize;
					int indexB = (secondCityIndex-j + genotypeSize) % genotypeSize; 

					Object temp = clonedIndividual.genotype.get(indexA);//store temp
					clonedIndividual.genotype.set(indexA, clonedIndividual.genotype.get(indexB));
					clonedIndividual.genotype.set(indexB, temp);
					if (Math.abs(indexA-indexB) <= 1 || Math.abs(indexA-indexB) >= maxIndex) break;
					
					j++;
				}

				firstCity = secondCity;
				// The new index for secondCity should be right next to the first city due to the swap
				firstCityIndex = (firstCityIndex+1) % genotypeSize;
			}
			
			if(Config.getInstance().calculateFitness(clonedIndividual)>=Config.getInstance().calculateFitness(originalIndividual)){//compare fitness'
				p.population.set(individualIndex, clonedIndividual);
			}
		}
		
		return p;
	}
	
	private int getIndexOfElement(Object element, List<Object> list){
		int index;
		for (index = 0; index < list.size(); index++){
			if (list.get(index).equals(element)){
				break;
			}
		}
		return index;
	}
}
