package org.pathvisio.merger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import org.bridgedb.AttributeMapper;
import org.bridgedb.BridgeDb;
import org.bridgedb.DataSource;
import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;
import org.bridgedb.bio.DataSourceTxt;
import org.pathvisio.core.model.ConverterException;
import org.pathvisio.core.model.MGroup;
import org.pathvisio.core.model.MLine;
import org.pathvisio.core.model.ObjectType;
import org.pathvisio.core.model.Pathway;
import org.pathvisio.core.model.PathwayElement;
import org.pathvisio.core.model.PathwayElement.MAnchor;
import org.pathvisio.merger.graph.Edge;
import org.pathvisio.merger.graph.Graph;
import org.pathvisio.merger.graph.Node;
import org.pathvisio.merger.graph.XGMMLWriter;

/**
 * 
 * @author mkutmon
 *
 */
public class PathwayMerger {

	private File directory;
	private File bridgedbGene;
	private File bridgedbMetabolites;
	private AttributeMapper geneAttr;
	private IDMapper geneMapper;
	private IDMapper metMapper;
	
	private static String PATHWAY_DIR = "pathway.dir";
	private static String OUTPUT_FILE = "output.file";
	private static String GENE_BRIDGEDB = "gene.bridgedb";
	private static String METABOLITE_BRIDGEDB = "metabolite.bridgedb";
	private static String LOG_FILE = "log.file";
	
	private static BufferedWriter log;
	
	public static void main (String [] args) throws Exception {
		if(args.length == 1) {
			File propsFile = new File(args[0]);
			if(propsFile.exists()) {
				props = new Properties();
				props.load(new FileReader(new File(args[0])));
				
				logFile = new File(props.getProperty(LOG_FILE));
				log = new BufferedWriter(new FileWriter(logFile));
				
				PathwayMerger gen = new PathwayMerger();
				gen.init();
				Graph graph = gen.createNetwork();
				
				File output = new File(props.getProperty(OUTPUT_FILE));
				XGMMLWriter.write(graph, new PrintWriter(output));
				log.close();
			}
		}
	}
	

	private static Properties props; 
	private static File logFile;

	public PathwayMerger() throws ClassNotFoundException {
		directory = new File(props.getProperty(PATHWAY_DIR));
		bridgedbGene = new File(props.getProperty(GENE_BRIDGEDB));
		bridgedbMetabolites = new File(props.getProperty(METABOLITE_BRIDGEDB));

		
		nodes = new HashMap<String, Node>();
		map = new HashMap<PathwayElement, Node>();
		edges = new HashMap<String, Edge>();
		graph = new Graph();
	}
	
	private List<Pathway> pathways;
	private Graph graph;
	private Map<String, Node> nodes;
	private Map<PathwayElement, Node> map;
	private Map<String, Edge> edges;
	
	public Graph createNetwork() throws IDMapperException, ConverterException, IOException {
		pathways = readPathways(directory);
		log.write("Parsing pathways from " + directory.getAbsolutePath() + "\n... containing " + pathways.size() + " pathways.");
		for(Pathway pathway : pathways) {
			log.write("\n> Parse pathway " + pathway.getMappInfo().getMapInfoName() + " with " + pathway.getDataObjects().size() + " pathway elements.\n");
			int genes = 0;
			int metabolites = 0;
			int pathways = 0;
			int groups = 0;
			int edges = 0;
			for(PathwayElement e : pathway.getDataObjects()) {
				if(e.getXref() != null && !e.getXref().getId().equals("") && e.getXref().getDataSource() != null) {
					if(e.getDataNodeType().equals("GeneProduct") || e.getDataNodeType().equals("Protein")) {
						// all genes and proteins are mapped to Ensembl
						genes++;
						createNode(e, geneMapper, "En", pathway);
					} else if(e.getDataNodeType().equals("Metabolite")) {
						// all metabolites are mapped to HMDB
						metabolites++;
						createNode(e, metMapper, "Ch", pathway);
					} else if(e.getDataNodeType().equals("Pathway")) {
						pathways++;
						createNode(e, null, "Wp", pathway);
					}
				} else {
					
				}
			}
			for(PathwayElement e : pathway.getDataObjects()) {
				if(e.getObjectType().equals(ObjectType.GROUP)) {
					groups++;
					createGroup(e, pathway);
				}
			}
			for(PathwayElement e : pathway.getDataObjects()) {
				if(e.getObjectType().equals(ObjectType.LINE)) {
					edges++;
					createEdge(e, pathway);
				}
			}
			log.write("\tGene count: " + genes + "\n");
			log.write("\tMetabolite count: " + metabolites + "\n");
			log.write("\tPathway count: " + pathways + "\n");
			log.write("\tGroup count: " + groups + "\n");
			log.write("\tEdge count: " + edges + "\n");
			log.write("\n\n");
		}
		
		for(Pathway p : pathways) {
			graph.appendAttribute(pathways.indexOf(p) + " Pathway", p.getMappInfo().getMapInfoName());
		}
		
		System.out.println("Conversion finished with " + nodes.size() + " nodes and " + edges.size() + " edges.");
		log.write("\n\nConversion finished with " + nodes.size() + " nodes and " + edges.size() + " edges.");
		return graph;
	}
	
	
	
	private void createGroup(PathwayElement e, Pathway pathway) {
		MGroup group = (MGroup)e;
		List<PathwayElement> list = new ArrayList<PathwayElement>();
		for(PathwayElement groupElement : group.getGroupElements()) {
			if(map.containsKey(groupElement)) {
				list.add(groupElement);
			}
		}
		if(list.size() >= 2) {
			String id = pathways.indexOf(pathway) + "." + e.getGroupId();
			Node groupNode = graph.addNode(id);
			groupNode.appendAttribute("Type", "Group");
			groupNode.appendAttribute("pathways", pathways.indexOf(pathway) +"");
			map.put(e, groupNode);
			nodes.put(id, groupNode);
			
			for(PathwayElement element : list) {
				String edgeId = map.get(element).getId() + " - " + id;
				if(!edges.containsKey(edgeId)) {
					Edge edge = graph.addEdge(edgeId, map.get(element), groupNode);
					edge.appendAttribute("Type", "Group");
					edge.appendAttribute("pathways", pathways.indexOf(pathway) +"");
					edges.put(edgeId, edge);
				} else {
					Edge edge = graph.getEdge(edgeId);
					String attr = (String) edge.getAttribute("pathways");
					if(!attr.contains(pathways.indexOf(pathway)+"")) {
						edge.setAttribute("pathways", attr + " | " + pathways.indexOf(pathway)+"");
					}
				}
			}
		}
	}

	private void createEdge(PathwayElement e, Pathway pathway) {
		MLine line = (MLine)e;
		PathwayElement start = pathway.getElementById(line.getStartGraphRef());
		PathwayElement end = pathway.getElementById(line.getEndGraphRef());
		if(line.getMAnchors().size() == 0) {
			if(start != null && end != null) {
				if(map.containsKey(start) && map.containsKey(end)) {
					String id = map.get(start).getId() + " - " + map.get(end).getId();
					if(!edges.containsKey(id)) {
						Edge edge = graph.addEdge(id, map.get(start), map.get(end));
						edge.appendAttribute("Type", line.getStartLineType().getName());
						edge.appendAttribute("pathways", pathways.indexOf(pathway) +"");
						edges.put(id, edge);
					} else {
						Edge edge = graph.getEdge(id);
						String attr = (String) edge.getAttribute("pathways");
						if(!attr.contains(pathways.indexOf(pathway)+"")) {
							edge.setAttribute("pathways", attr + " | " + pathways.indexOf(pathway)+"");
						}
					}
				}
			}
		} else {
			List<PathwayElement> list = new ArrayList<PathwayElement>();
			if(start != null && map.containsKey(start)) list.add(start);
			if(end != null && map.containsKey(end)) list.add(end);
			for(MAnchor anchor : line.getMAnchors()) {
				for(PathwayElement l : pathway.getDataObjects()) {
					if(l.getObjectType().equals(ObjectType.LINE)) {
						if(l.getStartGraphRef() != null && l.getEndGraphRef() != null) {
							if(l.getStartGraphRef().equals(anchor.getGraphId())) {
								if(map.containsKey(pathway.getElementById(l.getEndGraphRef()))) {
									list.add(pathway.getElementById(l.getEndGraphRef()));
								}
							} else if(l.getEndGraphRef().equals(anchor.getGraphId())) {
								if(map.containsKey(pathway.getElementById(l.getStartGraphRef()))) {
									list.add(pathway.getElementById(l.getStartGraphRef()));
								}
							}
						}
					}
				}
			}
			if(list.size() == 2) {
				String id = map.get(list.get(0)).getId() + " - " + map.get(list.get(1)).getId();
				if(!edges.containsKey(id)) {
					Edge edge = graph.addEdge(id, map.get(list.get(0)), map.get(list.get(1)));
					edge.appendAttribute("Type", "Anchor");
					edge.appendAttribute("pathways", pathways.indexOf(pathway) +"");
					edges.put(id, edge);
				} else {
					Edge edge = graph.getEdge(id);
					String attr = (String) edge.getAttribute("pathways");
					if(!attr.contains(pathways.indexOf(pathway)+"")) {
						edge.setAttribute("pathways", attr + " | " + pathways.indexOf(pathway)+"");
					}
				}
			} else if (list.size() > 2) {
				String id = pathways.indexOf(pathway) + "." + line.getMAnchors().get(0).getGraphId();
				Node anchorNode = graph.addNode(id);
				anchorNode.appendAttribute("Type", "Anchor");
				anchorNode.appendAttribute("pathways", pathways.indexOf(pathway) +"");
				map.put(line, anchorNode);
				nodes.put(id, anchorNode);
				
				for(PathwayElement element : list) {
					String edgeId = map.get(element).getId() + " - " + id;
					Edge edge = graph.addEdge(edgeId, map.get(element), anchorNode);
					edge.appendAttribute("Type", "Anchor");
					edge.appendAttribute("pathways", pathways.indexOf(pathway) +"");
					edges.put(edgeId, edge);
				}
			}
		}
	}
	
	private void createNode(PathwayElement e, IDMapper mapper, String systemCode, Pathway pathway) throws IDMapperException {
		if(mapper == null) {
			if(e.getXref() == null) {
				Node node = graph.addNode(e.getGraphId());
				node.appendAttribute("Label", e.getTextLabel());
				node.appendAttribute("pathways", pathways.indexOf(pathway) +"");
				node.appendAttribute("pathwayCount", 1+"");
				node.appendAttribute("Type", e.getDataNodeType());
				map.put(e, node);
			} else {
				if(!nodes.containsKey(e.getXref().getId())) {
					Node node = graph.addNode(e.getXref().getId());
					node.appendAttribute("GeneId", e.getXref().getId());
					node.appendAttribute("UnifiedId", "");
					node.appendAttribute("Label", e.getTextLabel());
					node.appendAttribute("pathways", pathways.indexOf(pathway) +"");
					node.appendAttribute("pathwayCount", 1+"");
					node.appendAttribute("Type", e.getDataNodeType());
					map.put(e, node);
					nodes.put(e.getXref().getId(), node);
				} else {
					Node node = nodes.get(e.getXref().getId());
					String attr = (String) node.getAttribute("pathways");
					int attrCount = Integer.parseInt((String)node.getAttribute("pathwayCount"));
					if(!attr.contains(pathways.indexOf(pathway)+"")) {
						node.setAttribute("pathwayCount", (attrCount+1)+"");
						node.setAttribute("pathways", attr + " | " + pathways.indexOf(pathway)+"");
					}
					map.put(e, node);
				}
			}
		} else {
			Xref unifiedId = null;
			String syscode;
			if(e.getXref().getDataSource().getFullName().equals("Uniprot/TrEMBL")) {
				syscode = "S";
			} else {
				syscode = e.getXref().getDataSource().getSystemCode();
			}
			if(syscode != null) {
				if(syscode.equals(systemCode)) {
					unifiedId = e.getXref();
				} else {
					Set<Xref> res = mapper.mapID(e.getXref(), DataSource.getExistingBySystemCode(systemCode));
					if(res.size() != 0) {
						unifiedId = res.iterator().next();
					}
				}
			}
			
			if(unifiedId != null) {
				if(!nodes.containsKey(unifiedId.getId())) {
					Node node = graph.addNode(unifiedId.getId());
					node.appendAttribute("GeneId", e.getXref().getId());
					node.appendAttribute("UnifiedId", unifiedId.getId());
					node.appendAttribute("Label", getLabel(e, unifiedId));
					node.appendAttribute("pathways", pathways.indexOf(pathway) +"");
					node.appendAttribute("pathwayCount", 1+"");
					node.appendAttribute("Type", e.getDataNodeType());
					map.put(e, node);
					nodes.put(unifiedId.getId(), node);
				} else {
					Node node = nodes.get(unifiedId.getId());
					String attr = (String) node.getAttribute("pathways");
					int attrCount = Integer.parseInt((String)node.getAttribute("pathwayCount"));
					if(!attr.contains(pathways.indexOf(pathway)+"")) {
						node.setAttribute("pathwayCount", (attrCount+1)+"");
						node.setAttribute("pathways", attr + " | " + pathways.indexOf(pathway)+"");
					}
					map.put(e, node);
				}
			} else {
				String id = e.getXref().getId();
				if(!nodes.containsKey(id)) {
					Node node = graph.addNode(id);
					node.appendAttribute("GeneId", e.getXref().getId());
					node.appendAttribute("Label", getLabel(e, e.getXref()));
					node.appendAttribute("pathways", pathways.indexOf(pathway) +"");
					node.appendAttribute("pathwayCount", 1+"");
					node.appendAttribute("Type", e.getDataNodeType());
					map.put(e, node);
					nodes.put(id, node);
				} else {
					Node node = nodes.get(id);
					String attr = (String) node.getAttribute("pathways");
					int attrCount = Integer.parseInt((String)node.getAttribute("pathwayCount"));
					if(!attr.contains(pathways.indexOf(pathway)+"")) {
						node.setAttribute("pathwayCount", (attrCount+1)+"");
						node.setAttribute("pathways", attr + " | " + pathways.indexOf(pathway)+"");
					}
					map.put(e, node);
				}
			}
		}
	}
	
	private String getLabel(PathwayElement element, Xref xref) throws IDMapperException {
		if(element.getDataNodeType().equals("Metabolite")) {
			return element.getTextLabel();
		} else {
			Map<String, Set<String>> attrMap = geneAttr.getAttributes(xref);
			Set<String> symbol = attrMap.get("Symbol");
			if(symbol != null && symbol.size() > 0) {
				return symbol.iterator().next();
			} else {
				return element.getTextLabel();
			}
		}
	}
	
	private List<Pathway> readPathways(File directory) throws ConverterException {
		List<Pathway> list = new ArrayList<Pathway>();
		
		for(File file : directory.listFiles()) {
			if(file.getName().endsWith(".gpml")) {
				Pathway pathway = new Pathway();
				pathway.readFromXml(file, true);
				list.add(pathway);
			}
		}
		
		return list;
	}
	
	private void init() throws IDMapperException, ClassNotFoundException {
		DataSourceTxt.init();
		Class.forName("org.bridgedb.rdb.IDMapperRdb");
		
		geneMapper = BridgeDb.connect("idmapper-pgdb:" + bridgedbGene.getAbsolutePath());
		metMapper = BridgeDb.connect("idmapper-pgdb:" + bridgedbMetabolites.getAbsolutePath());
		geneAttr = (AttributeMapper) geneMapper;
	}
}
