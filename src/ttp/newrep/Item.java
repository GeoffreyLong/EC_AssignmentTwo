package ttp.newrep;

public class Item {
	public int weight;
	public int profit;
	public boolean isSelected;
	
	public Item(int profit, int weight){
		this.profit = profit;
		this.weight = weight;
		this.isSelected = false;
	}
}
