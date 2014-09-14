package ttp.newrep;

public class Individual {
	public City[] tour;
	public City startingCity;
	
	public Individual(double[][] nodes, int[][] items, int[] tour, int itemsPerCity){
		this.tour = new City[tour.length-2];
		startingCity = new City(nodes[tour[0]][0],nodes[tour[0]][1],nodes[tour[0]][2]);
		for (int index = 1; index < tour.length-1; index ++){
			int tourIndex = tour[index];
			City city = new City(nodes[tourIndex][0], nodes[tourIndex][1], nodes[tourIndex][2]);
			for (int i = 0; i < itemsPerCity; i++){
				int itemIndex = (tour.length-2) * i + tourIndex - 1;
				Item item = new Item(items[itemIndex][1], items[itemIndex][2]);
				city.items.add(item);
			}
			this.tour[index-1] = city;
		}
	}
}
