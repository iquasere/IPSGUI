package interProParser;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Parser {
	
	public static String jobID = "iprscan5-R20160424-232342-0215-12126963-pg";
	public static String directory = "C:\\Users\\João Sequeira\\Documents\\Escola\\Projeto Bioinformática BSIstemas\\interproscan API\\";

	public static ArrayList<String> readFile(String jobID, String extension) throws IOException{
		
		Scanner s = new Scanner(new File(directory+jobID+extension));
		ArrayList<String> list = new ArrayList<String>();
		while (s.hasNextLine()){
		    list.add(s.nextLine());
		}
		s.close();
		
		return list;
	}

	public static String get_protein_sequence() throws IOException{
		
		String protein_sequence = "";
		
		ArrayList<String> file = readFile(jobID,".gff.txt");
		
		for (int i = 0; i < file.size(); i ++) {
			
			if (file.get(i).startsWith(">EMBOSS_001")) {
				
				i += 1;
				
				while (!(file.get(i).startsWith(">match$"))){
				
					protein_sequence += file.get(i);
					
					i += 1;
					
				}
			}
		}
		
		return protein_sequence;	
		
	}
	
	public static List<HashMap<String, String>> get_gff_information() throws IOException {
		
		ArrayList<String> file = readFile(jobID,".gff.txt");
		
		List<HashMap<String, String>> domains = new ArrayList<HashMap<String, String>>();

		for (int i = 0; i < file.size(); i++) {
			
			if (file.get(i).startsWith("EMBOSS_001\t")) {
				
				if (!(file.get(i).startsWith("EMBOSS_001\t."))) {
					
					String[] parts = file.get(i).split("\t");
					HashMap<String,String> domain = new HashMap<String,String>();
					
					domain.put("Source", parts[1]);      			//database from where the result is
					domain.put("Type", parts[2]);					//type of match
					domain.put("Start_position", parts[3]);
					domain.put("End_position", parts[4]);
					domain.put("E-value", parts[5]);
					domain.put("Strand", parts[6]);
					domain.put("Phase", parts[7]);					//phase indicates where the feature begins with reference to the reading frame
					
					Pattern GOPattern = Pattern.compile("Ontology_term=(.+?);");
					Matcher GOMatcher = GOPattern.matcher(parts[8]);
					if (GOMatcher.find()){domain.put("Ontology_term", GOMatcher.group(1));}
					
					Pattern IDPattern = Pattern.compile("ID=match(.+)_[0-9]+_[0-9]+;");
					Matcher IDMatcher = IDPattern.matcher(parts[8]);
					if(IDMatcher.find()){domain.put("Match_ID", IDMatcher.group(1));}
					
					Pattern SDPattern = Pattern.compile("signature_desc=(.+?);");
					Matcher SDMatcher = SDPattern.matcher(parts[8]);
					if(SDMatcher.find()){domain.put("Description", SDMatcher.group(1));}
					
					Pattern NamePattern = Pattern.compile("Name=(.+?);");
					Matcher NameMatcher = NamePattern.matcher(parts[8]);
					if(NameMatcher.find()){domain.put("Name", NameMatcher.group(1));}
					
					Pattern StatusPattern = Pattern.compile("status=(.+?);");
					Matcher StatusMatcher = StatusPattern.matcher(parts[8]);
					if(StatusMatcher.find()){domain.put("Status", StatusMatcher.group(1));}
					
					Pattern xrefPattern = Pattern.compile("Dbxref=(.+)");
					Matcher xrefMatcher = xrefPattern.matcher(parts[8]);
					if(xrefMatcher.find()){domain.put("Dbxref", xrefMatcher.group(1));}
					
					domains.add(domain);
					
				}
				
			}
				
			if (file.get(i).startsWith(">match")) {
							
				Pattern SeqPattern = Pattern.compile(">match(.+?)_[0-9]+_[0-9]+");
				Matcher SeqMatcher = SeqPattern.matcher(file.get(i));
				SeqMatcher.matches();
				
				String sequence	= "";
				
				i++;
				
				while (!(file.get(i).startsWith(">match")) && i < file.size()-1) {
					
					sequence += file.get(i);		
					
					i ++;
								
				}
				
				for (HashMap<String,String>domain:domains){
					
					if(domain.get("Match_ID").equals(SeqMatcher.group(1))){
						
						domain.put("Sequence",sequence);
					}
				}
				i--;
			}
		}
		
		return domains;	
	}
	
	public static List<HashMap<String, Object>> get_xml_information() throws IOException {
		
		ArrayList<String> file = readFile(jobID,".xml.xml");
		
		List<HashMap<String, Object>> domains = new ArrayList<HashMap<String, Object>>();
		
		Pattern matchPattern = Pattern.compile("<([^/]+)-match");
		Pattern descPattern = Pattern.compile("desc=\"(.+?)\"");
		Pattern typePattern = Pattern.compile("type=\"(.+?)\"");
		Pattern evaluePattern = Pattern.compile("evalue=\"(.+?)\"");
		Pattern scorePattern = Pattern.compile("score=\"(.+?)\"");
		Pattern familyNamePattern = Pattern.compile("familyName=\"(.+?)\"");
		Pattern acPattern = Pattern.compile("ac=\"(.+?)\"");
		Pattern namePattern = Pattern.compile("name=\"(.+?)\"");
		Pattern xrefPattern = Pattern.compile("\"(.+?)\"");
		Pattern modelPattern = Pattern.compile("\"(.+?)\"");
		Pattern libraryPattern = Pattern.compile("library=\"(.+?)\"");
		Pattern versionPattern = Pattern.compile("version=\"(.+?)\"");
		Pattern hmmStartPattern = Pattern.compile("hmm-start=\"(.+?)\"");
		Pattern hmmEndPattern = Pattern.compile("hmm-end=\"(.+?)\"");
		Pattern lengthPattern = Pattern.compile("hmm-length=\"(.+?)\"");
		Pattern startPattern = Pattern.compile("start=\"(.+?)\"");
		Pattern endPattern = Pattern.compile("end=\"(.+?)\"");
		
		for (int i = 6; i < file.size(); i++) {
			
			Matcher matchMatcher = matchPattern.matcher(file.get(i));
			
			if (matchMatcher.find()){
				
				HashMap<String,Object> domain = new HashMap<String,Object>();
				
				Matcher evalueMatcher = evaluePattern.matcher(file.get(i));
				if (evalueMatcher.find()){domain.put("E-value",Float.parseFloat(evalueMatcher.group(1)));}
				
				Matcher scoreMatcher = scorePattern.matcher(file.get(i));
				if (scoreMatcher.find()){domain.put("Score",Float.parseFloat(scoreMatcher.group(1)));}
				
				Matcher familyNameMatcher = familyNamePattern.matcher(file.get(i));
				if (familyNameMatcher.find()){domain.put("Family_Name",familyNameMatcher.group(1));}
				
				i++;
		
				Matcher acMatcher = acPattern.matcher(file.get(i));
				if (acMatcher.find()){domain.put("Signature_accession_number",acMatcher.group(1));}
				
				Matcher nameMatcher = namePattern.matcher(file.get(i));
				if (nameMatcher.find()){domain.put("Name",nameMatcher.group(1));}
				
				i++;
				
				if (file.get(i).trim().startsWith("<entry")) {
					
					Matcher entryAcMatcher = acPattern.matcher(file.get(i));
					if (entryAcMatcher.find()) {domain.put("Entry_accession_number",entryAcMatcher.group(1));}
					
					Matcher descMatcher = descPattern.matcher(file.get(i));
					if (descMatcher.find()) {domain.put("Entry_description", descMatcher.group(1));}
					
					Matcher entryNameMatcher = namePattern.matcher(file.get(i));
					if (entryNameMatcher.find()) {domain.put("Entry_name", entryNameMatcher.group(1));}
					
					Matcher typeMatcher = typePattern.matcher(file.get(i));
					if (typeMatcher.find()) {domain.put("Entry_type", entryNameMatcher.group(1));}
					
					i++;
					
					List<List<String>> allXref = new ArrayList<List<String>>();
					
					//entries with Xrefs end in </entry>, but those that lack xRefs end at the start of the models
					while(!(file.get(i).trim().startsWith("</entry>")) && !(file.get(i).trim().startsWith("<models>")) && i < file.size()) {
						 Matcher xrefMatcher = xrefPattern.matcher(file.get(i));
						 
						 List<String> xref = new ArrayList<String>();
						 
						 while (xrefMatcher.find()) {

							 xref.add(xrefMatcher.group().replace("\"", ""));
						 }
						 
						 allXref.add(xref);
							
						 i++;
					}
					
					if (allXref.size() > 0) {
						
						domain.put("Xrefs", allXref);
						
						i++;
					}
				}
				
				if (file.get(i).trim().startsWith("<models>")) {
					
					i++;
				
					List<List<String>> allModels = new ArrayList<List<String>>();
					
					while(!(file.get(i).trim().startsWith("</models>")) && i < file.size()) {
						 Matcher modelMatcher = modelPattern.matcher(file.get(i));
						 
						 List<String> Model = new ArrayList<String>();
						 
						 while (modelMatcher.find()) {

							 Model.add(modelMatcher.group().replace("\"", ""));
						 }
						 
						 allModels.add(Model);
							
						 i++;
					}
					
					domain.put("Models", allModels);
					
					i++;
					
				}
				
				Matcher libraryMatcher = libraryPattern.matcher(file.get(i));
				if (libraryMatcher.find()){domain.put("Library",libraryMatcher.group(1));}
				
				Matcher versionMatcher = versionPattern.matcher(file.get(i));
				if (versionMatcher.find()){domain.put("Library_version",versionMatcher.group(1));}
					
				i += 3;
				
				List<HashMap<String,Float>> allLocations = new ArrayList<HashMap<String,Float>>();

				while (!(file.get(i).trim().startsWith("</locations>"))){
					
					HashMap<String,Float> location = new HashMap<String,Float>();
					
					Matcher locationScoreMatcher = scorePattern.matcher(file.get(i));
					if (locationScoreMatcher.find()){location.put("Score",Float.parseFloat((locationScoreMatcher.group(1))));}
					
					Matcher locationEvalueMatcher = evaluePattern.matcher(file.get(i));
					if (locationEvalueMatcher.find()){location.put("E-value",Float.parseFloat((locationEvalueMatcher.group(1))));}
					
					Matcher hmmStartMatcher = hmmStartPattern.matcher(file.get(i));
					if (hmmStartMatcher.find()){location.put("HMM-Start",Float.parseFloat((hmmStartMatcher.group(1))));}
					
					Matcher hmmEndMatcher = hmmEndPattern.matcher(file.get(i));
					if (hmmEndMatcher.find()){location.put("HMM-End",Float.parseFloat(hmmEndMatcher.group(1)));}
					
					Matcher lengthMatcher = lengthPattern.matcher(file.get(i));
					if (lengthMatcher.find()){location.put("HMM-Length",Float.parseFloat(lengthMatcher.group(1)));}
				
					Matcher startMatcher = startPattern.matcher(file.get(i));
					if (startMatcher.find()){location.put("Start",Float.parseFloat(startMatcher.group(1)));}
					
					Matcher endMatcher = endPattern.matcher(file.get(i));
					if (endMatcher.find()){location.put("End",Float.parseFloat(endMatcher.group(1)));}
					
					if(location.size()>0) {allLocations.add(location);}
					
					i++;
				
				}
				
				domain.put("Locations", allLocations);
					
				domains.add(domain);
			}
		}	
		
		return domains;	
	}
	
	public static void main(String[] args) throws IOException{
		
		for (HashMap<String, Object>domain:get_xml_information()){
			
//			List<HashMap<String,Float>> Locations = new ArrayList<HashMap<String,Float>> (domain.get("Locations"));
		
			
			System.out.println(domain);
		}
	}
}
