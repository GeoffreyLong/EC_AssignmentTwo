package ttp.newrep;

public class Item {
	public int weight;
	public int profit;
	public boolean isSelected;
	public int itemId;
	
	public Item(int profit, int weight, int itemId){
		this.profit = profit;
		this.weight = weight;
		this.itemId = itemId;
		this.isSelected = false;
	}
}
