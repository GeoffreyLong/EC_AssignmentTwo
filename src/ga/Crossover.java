package ga;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;


public class Crossover {
	private static final boolean verbose = false;
	private static final Random rand = new Random(System.currentTimeMillis());
	private CrossoverType crossoverType;
	private Selection pSelect;
	private double[] crossoverTypeChance;
	
	public enum CrossoverType{
		ORDER, PMX, CYCLE, EDGE;
	}
	
	public Crossover(double[] crossoverTypeChance){
		this.crossoverTypeChance = crossoverTypeChance;
		pSelect = new Selection(Config.getInstance().getParentSelectionType());
	}
	public Population cross(Population p){
	
		Config conf = Config.getInstance();
		
		// Select mating pool
		Population matingPool = pSelect.select(p);
		
		// Apply crossover to pairs
		Population offspring = new Population();
		
		double random = rand.nextDouble();
		
		crossoverType = null;
		for (int i = 0; i < CrossoverType.values().length; i++){
			if (random < crossoverTypeChance[i]){
				crossoverType = CrossoverType.values()[i];
				break;
			}
		}

		if (crossoverType == null){
			return p;
		}
		
		for (int i = 0; i < matingPool.size()/2; i++) {
			Individual p1 = matingPool.population.get(i*2);
			Individual p2 = matingPool.population.get(i*2+1);			
			
			if (rand.nextDouble() < conf.crossoverChance) {
				offspring.population.add(cross(p1,p2));
				offspring.population.add(cross(p2, p1));
			} else {
				offspring.population.add(p1);
				offspring.population.add(p2);
			}			
		}
		if ( (matingPool.size() % 2) == 1 ) {
			Individual ind = matingPool.population.get(matingPool.size()-1);
			offspring.population.add(ind);
		}
			
		
		return offspring;
		
	}
	
	public Individual cross(Individual a, Individual b) {
		switch(crossoverType) {
			case ORDER: 
				return this.orderCross(a, b);
			case PMX:
				return this.pmxCross(a, b);				
			case CYCLE:
				return this.cycleCross(a, b);
			case EDGE:
				return this.edgeRecombination(a, b);
			default: throw new IllegalArgumentException("Crossover type not specified");				
		}	
	}
	
	public Individual orderCross(Individual a, Individual b) {
		int numChromosomes = a.genotype.size();
		if (verbose) System.out.println(numChromosomes);
		
		
		List<Object> offspringGenotype = new ArrayList<Object>(Collections.nCopies(numChromosomes, -1));
		
		// Generate random cut interval
		int startCut = rand.nextInt(numChromosomes);
		int endCut = rand.nextInt(numChromosomes);
		/*
		if (startCut > endCut) {
			int tmp = startCut;
			startCut = endCut;
			endCut = tmp;
		}
		*/
		if (verbose) {
			System.out.printf("startCut:%d, endCut:%d\n",startCut,endCut);
			System.out.println("before copy");
			System.out.println(offspringGenotype);
		}
		// Copy the cut sections to the new individuals, wrapping over at end
		int numCopied = 1;
		offspringGenotype.set(endCut,a.genotype.get(endCut));
		for (int i = startCut; i != endCut; i = (i + 1) % numChromosomes) {
			offspringGenotype.set(i,a.genotype.get(i));
			numCopied++;
		}
	
		if (verbose) {
			System.out.println("after copy");
			System.out.println(offspringGenotype);
		}
		// Copy the remaining chromosomes
		// TODO: contains in individual probably won't work, need to implement equals etc\
		int chromosomesRemaining = numChromosomes - numCopied;
		int setIndex = (endCut + 2 > numChromosomes) ? 0 : endCut + 1;
		int getIndex = setIndex;
		while (chromosomesRemaining > 0) {
			Object chromosome = b.genotype.get(getIndex);
			if (!offspringGenotype.contains(chromosome)) {
				offspringGenotype.set(setIndex, chromosome);
				setIndex = (setIndex + 2 > numChromosomes) ? 0 : setIndex + 1;
				chromosomesRemaining--;
			}
			getIndex = (getIndex + 2 > numChromosomes) ? 0 : getIndex + 1;
		}
		
		return new Individual(offspringGenotype);
	}
	
	public Individual pmxCross(Individual a, Individual b) {
		int numChromosomes = a.genotype.size();
		
		List<Object> offspringGenotype = new ArrayList<Object>(Collections.nCopies(numChromosomes, -1));
		
		// Generate random cut interval
		int startCut = rand.nextInt(numChromosomes);
		int endCut = rand.nextInt(numChromosomes);
		
		if (verbose) {
			System.out.printf("startCut:%d, endCut:%d\n",startCut,endCut);
			System.out.println("before copy");
			System.out.println(offspringGenotype);
		}
		// Copy the cut sections to the new individuals
		List<Integer> copyIndex = new ArrayList<Integer>();
		for (int i = startCut; i != endCut; i = (i + 1) % numChromosomes) {
			copyIndex.add(i);
		}
		copyIndex.add(endCut);
		for (Integer i : copyIndex) {
			offspringGenotype.set(i,a.genotype.get(i));
		}
		
		if (verbose) {
			System.out.println("after copy");
			System.out.println(offspringGenotype);
		}
		
		// Construct element index lookup table
		Map<Object, Integer> bLookup = new HashMap<Object, Integer>();
		{int index = 0;
		for (Object chromosome : b.genotype) {
			bLookup.put(chromosome, index);
			index++;
		}}
		
		for (Integer i : copyIndex) {
			Object el1 = b.genotype.get(i);
			if (!offspringGenotype.contains(el1)) {
				boolean placeFound = false;
				int index = i;
				Object el2 = offspringGenotype.get(index);
				while (!placeFound) {					
					// find index of el2 in b
					index = bLookup.get(el2);
					el2 = offspringGenotype.get(index);
					placeFound = (el2 instanceof Integer) ? ((int)el2 == -1) : false;
				}
				offspringGenotype.set(index,el1);
			}
		}
		if (verbose) {
			System.out.println("after pmx");
			System.out.println(offspringGenotype);
		}
		
		for (int i = 0; i < offspringGenotype.size(); i++) {
			Object el = offspringGenotype.get(i);
			boolean placeFilled = (el instanceof Integer) ? ((int)el != -1) : true;
			if (verbose) System.out.printf("index: %d, placeFilled:%s, integer:%s\n", i, placeFilled, el instanceof Integer);
			if (!placeFilled) {
				offspringGenotype.set(i, b.genotype.get(i));
			}
		}
		
		if (verbose) {
			System.out.println("after direct mapping");
			System.out.println(offspringGenotype);
		}
		return new Individual(offspringGenotype);
		
	}
	
	public Individual cycleCross(Individual a, Individual b) {
		// Construct element index lookup table
		Map<Object, Integer> aLookup = new HashMap<Object, Integer>();
		{int index = 0;
		for (Object chromosome : a.genotype) {
			aLookup.put(chromosome, index);
			index++;
		}}
		
		// Identify cycles
		List<List<Integer>> cycleList = new ArrayList<List<Integer>>();
		Set<Integer> indexDone = new HashSet<Integer>();
	
		int startIndex = 0;
		while (indexDone.size() < a.genotype.size()) { // Loop until all index's are in a cycle
			if (indexDone.contains(startIndex)) { // find next index
				startIndex++;
				continue;
			}
			if (verbose) System.out.println(startIndex);
			// loop init
			List<Integer> cycle = new ArrayList<Integer>();
			Object cycleStart = a.genotype.get(startIndex);
			int index = startIndex;
					
			do { // follow a cycle, store indexs
				cycle.add(index);
				indexDone.add(index);
				index = aLookup.get(b.genotype.get(index));
				
			} while (!a.genotype.get(index).equals(cycleStart));
			cycleList.add(cycle);
			
			startIndex++;
		}
		
		// Debug cycle indexs
		if (verbose) {
			System.out.println("#####################");
			System.out.println("AGenotype: " + a.genotype);
			System.out.println("BGenotype: " + b.genotype);
			System.out.println("Generated Cycles");
			for (List<Integer> cycle : cycleList) {
				System.out.println(cycle);
			}
			System.out.println("#####################");
		}
		
		
		// Add cycle elements to offspring. Odd cycles add odd elements, even cycles add even elements.
		int numChromosomes = a.genotype.size();
		List<Object> offspringGenotype = new ArrayList<Object>(Collections.nCopies(numChromosomes, -1));
		boolean swap = true;
		for (List<Integer> cycle : cycleList) {
			for (int idx : cycle) {
				Object el = (swap) ? a.genotype.get(idx) : b.genotype.get(idx);
				offspringGenotype.set(idx, el);
			}
			swap = !swap;
		}
		
		
		
		return new Individual(offspringGenotype);
	}
	
	public Individual edgeRecombination(Individual a, Individual b) {
		List<Object> offspringGenotype = new ArrayList<Object>();
		
		Map<Object, List<Object>> edgeTable = new HashMap<Object, List<Object>>();
		for (Object el : a.genotype){
			edgeTable.put(el, new ArrayList<Object>(4));
		}
		
		// Populate edge table
		for (int i = 0; i < a.genotype.size(); i++) {
			int idxFwd = (i+1 >= a.genotype.size()) ? 0 : i+1;
			int idxBwd = (i-1 < 0) ? a.genotype.size() - 1 : i-1;
		
			edgeTable.get(a.genotype.get(i)).add(a.genotype.get(idxFwd));
			edgeTable.get(a.genotype.get(i)).add(a.genotype.get(idxBwd));
			edgeTable.get(b.genotype.get(i)).add(b.genotype.get(idxFwd));
			edgeTable.get(b.genotype.get(i)).add(b.genotype.get(idxBwd));
		}
		
		//debug print
		if (verbose) {
			System.out.println("Edge Table");
			for (Map.Entry<Object, List<Object>> entry : edgeTable.entrySet()) {
			    Object key = entry.getKey();
			    List<Object> value = entry.getValue();
			    System.out.println(key + ": " + value);
			}
		}
		
		Object startEdge = a.genotype.get(0);
		Set<Object> doneEdges = new HashSet<Object>();
		while (doneEdges.size() < a.genotype.size()) {
			if (edgeTable.get(startEdge).size() <= 0) {
				// TODO FIND NEW EDGE
				boolean found = false;
				for (Object edge : a.genotype) {
					if (edgeTable.get(edge).size() > 0) {						
						startEdge = edge;
						found = true;
						break;
					}
				}
				if (!found) {
					if (verbose) System.out.println("WHAT THIS SHOULDNT HAPPEN");
				}
			}
			// Count the double edges
			List<Object> currentList = edgeTable.get(startEdge);
			Map<Object,Integer> numOccurences = new HashMap<Object,Integer>();
			for (Object edge : currentList) {
				if (numOccurences.get(edge) == null) {					
					numOccurences.put(edge,1);
				} else {
					numOccurences.put(edge, 2); // MAX is 2
				}
			}
			// Get the double edges
			List<Object> doubleEdges = new ArrayList<Object>();
			for (Object edge : currentList) {
				if (numOccurences.get(edge) >= 2) {
					doubleEdges.add(edge);
				}
			}
			
			Object chosenEdge = null;
			if (doubleEdges.size() == 0) {
				// Compare based on list size
				// TODO  This seems broken... Index out of bound here
				chosenEdge = currentList.get(0);
				int listSize = edgeTable.get(chosenEdge).size();
				for (Object edge : currentList) {
					int newListSize = edgeTable.get(edge).size();
					if (verbose) System.out.println(newListSize);
					if (newListSize < listSize) {
						chosenEdge = edge;
						listSize = newListSize;
					}
				}				
			} else if (doubleEdges.size() >= 2) {
				// chose one randomly				
				chosenEdge = doubleEdges.get(rand.nextInt(2));
				
			} else {
				// only one double edge
				chosenEdge = doubleEdges.get(0);
			}
			
			// Remove chosen edge from all lists
			for (List<Object> edgeList : edgeTable.values()) {
				while (edgeList.remove(chosenEdge));
			}
			
			offspringGenotype.add(chosenEdge);
			doneEdges.add(chosenEdge);
			startEdge = chosenEdge;
			//debug print	
			if (verbose) {
				System.out.println("Current Solution: " + offspringGenotype);
				System.out.println("Edge Table");
				for (Map.Entry<Object, List<Object>> entry : edgeTable.entrySet()) {
				    Object key = entry.getKey();
				    List<Object> value = entry.getValue();
				    System.out.println(key + ": " + value);
				}
			}
		}
		
		return new Individual(offspringGenotype);
	}
}
