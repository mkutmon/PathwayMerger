package org.pathvisio.merger.graph;

/**
 * @author mkutmon
 * adapted from Thomas Kelder
 */
public class Node extends AttributeHolder {
	private String id;
	
	public Node(String id) {
		this.id = id;
	}
	
	public String getId() {
		return id;
	}
	
	public int hashCode() {
		return id.hashCode();
	}
}