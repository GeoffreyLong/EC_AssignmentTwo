package ttp.newrep;

import java.awt.geom.Point2D;
import java.util.LinkedList;
import java.util.List;

public class City {
	public List<Item> items = new LinkedList<Item>();
	public Point2D location;
	public int cityId;
	
	public City(double id, double x, double y){
		this.cityId = (int) id;
		this.location = new Point2D.Double(x,y);
	}
	
	public int getWeight(){
		int weight = 0;
		for (Item item : items){
			if(item.isSelected){
				weight += item.weight;
			}
		}
		return weight;
	}
	
	public int getProfit(){
		int profit = 0;
		for (Item item : items){
			if(item.isSelected){
				profit += item.profit;
			}
		}
		return profit;
	}
}
