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
		
		for (int i = 6; i < file.size(); i++) {
			
			Matcher matchMatcher = matchPattern.matcher(file.get(i));
			
			if (matchMatcher.find()){
				
				HashMap<String,Object> domain = new HashMap<String,Object>();
				
				Pattern evaluePattern = Pattern.compile("evalue=\"(.+?)\"");
				Matcher evalueMatcher = evaluePattern.matcher(file.get(i));
				if (evalueMatcher.find()){domain.put("E-value",evalueMatcher.group(1));}
				
				Pattern scorePattern = Pattern.compile("score=\"(.+?)\"");
				Matcher scoreMatcher = scorePattern.matcher(file.get(i));
				if (scoreMatcher.find()){domain.put("Score",scoreMatcher.group(1));}
				
				Pattern familyNamePattern = Pattern.compile("familyName=\"(.+?)\"");
				Matcher familyNameMatcher = familyNamePattern.matcher(file.get(i));
				if (familyNameMatcher.find()){domain.put("Family_Name",familyNameMatcher.group(1));}
				
				i++;
		
				Pattern acPattern = Pattern.compile("ac=\"(.+?)\"");
				Matcher acMatcher = acPattern.matcher(file.get(i));
				if (acMatcher.find()){domain.put("Signature_accession_number",acMatcher.group(1));}
				
				Pattern namePattern = Pattern.compile("name=\"(.+?)\"");
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
				}
				
				i++;
				
				int f = 1;
				
				while(f > 0 && i < file.size()) {
					
					Matcher modelAcMatcher = acPattern.matcher(file.get(i));
					if 
					
					if (file.get(i).trim().startsWith("</models>")) {f--;}
				}
				
				if (file.get(i).trim().startsWith("<signature-library-release")) {
					
					Pattern libraryPattern = Pattern.compile("library=\"(.+?)\"");
					Matcher libraryMatcher = libraryPattern.matcher(file.get(i));
					if (libraryMatcher.find()){domain.put("Library",libraryMatcher.group(1));}
					
					Pattern versionPattern = Pattern.compile("version=\"(.+?)\"");
					Matcher versionMatcher = versionPattern.matcher(file.get(i));
					if (versionMatcher.find()){domain.put("Library_version",versionMatcher.group(1));}
				}
					
				i++;
					
				}
				
				System.out.println(domain);
				
				domains.add(domain);
				
			}
			
			
		}				
		
		return domains;	
	}
	
	public static void main(String[] args) throws IOException{
		
		System.out.println(get_xml_information());
		
	}
}
