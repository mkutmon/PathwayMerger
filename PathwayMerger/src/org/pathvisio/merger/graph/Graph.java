package org.pathvisio.merger.graph;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * @author mkutmon
 * adapted from Thomas Kelder
 */
public class Graph extends AttributeHolder {
	private String title = "";
	
	private Map<String, Node> nodes = new HashMap<String, Node>();
	private Map<String, Edge> edges = new HashMap<String, Edge>();

	public Node addNode(String id) {
		Node n = nodes.get(id);
		if(n == null) { 
			n = new Node(id);
			nodes.put(id, n);
		}
		return n;
	}
	
	public Edge addEdge(String id, Node source, Node target) {
		Edge e = edges.get(id);
		if(e == null) {
			e = new Edge(id, source, target);
			edges.put(id, e);
		}
		return e;
	}
	
	// SETTERS & GETTERS
	
	public void setTitle(String title) { this.title = title; }
	public String getTitle() { return title; }
	
	public Node getNode(String id) { return nodes.get(id); }
	public Edge getEdge(String id) { return edges.get(id); }
	
	public Collection<Node> getNodes() { return nodes.values(); }
	public Collection<Edge> getEdges() { return edges.values(); }
}
