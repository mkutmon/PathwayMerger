package org.pathvisio.merger.graph;

/**
 * 
 * @author mkutmon
 * adapted from Thomas Kelder
 */
public class Edge extends AttributeHolder {
	private String id;
	private Node source;
	private Node target;
		
	public Edge(String id, Node source, Node target) {
		this.id = id;
		this.source = source;
		this.target = target;
	}
	
	public Node getSource() {
		return source;
	}
		
	public Node getTarget() {
		return target;
	}
		
	public String getId() {
		return id;
	}
}