import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class Read {
	
	public static String[] line;
	public static String id = "";
	
	public static Double AA = 0.0;
	public static Double AC = 0.0;
	public static Double AG = 0.0;
	public static Double AT = 0.0;
	public static Double CC = 0.0; 
	public static Double CG = 0.0;
	public static Double CT = 0.0;
	public static Double GG = 0.0;
	public static Double GT = 0.0;
	public static Double TT = 0.0;
			
	public static Double A = 0.0;
	public static Double C = 0.0;
	public static Double T = 0.0;
	public static Double G = 0.0;
	
	public static Double sum = 0.0;
	public static int minimum_depth = 5;
	public static double minimum_af = 0.00;
	
	private String CHR = "CHR";
	private String POS = "POS";
	private String REF = "REF";
	private String ALT = "ALT";
	private String GENOTYPE = "GENOTYPE";
	private String FORMAT_GT = "FORMAT_GT";
	private String FORMAT_DP = "FORMAT_DP";
	private String FORMAT_AD = "FORMAT_AD";
	private String FORMAT_GQ = "FORMAT_GQ";

	private String NAMES = "NAMES";
	private String string_split = "-";

	public HashMap<String, ArrayList<Double>> read_GTCounts(String gt_count, HashMap<String, String> reference_genotype) throws IOException, InterruptedException{
		
		/**HashMap that lists the allele frequencies at each position:**/
		HashMap<String, ArrayList<Double>> alleles = new HashMap<String, ArrayList<Double>>();
		
		/**Open genotype count file:**/
		FileReader file = null;
		try {
			file = new FileReader(gt_count);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		BufferedReader in = new BufferedReader(file);

		/**Read genotype count file:**/
		String row = null;
		while((row = in.readLine()) != null){
			
			/**Skip command lines:**/
			if(row.startsWith("#"))	continue;
			
			/**Split current line:**/
			line = row.split("\t");
			
			/**Initialize position in vector**/
			id = line[0]+"\t"+line[1];
			
			/**Initialize Reference Genotype:**/
			if(!reference_genotype.containsKey(id)) reference_genotype.put(id, line[2]+","+line[2]);
			
			/**Genotype Counts:**/
			AA = Double.valueOf(line[4]);
			AC = Double.valueOf(line[5]);
			AG = Double.valueOf(line[6]);
			AT = Double.valueOf(line[7]);
			CC = Double.valueOf(line[8]);
			CG = Double.valueOf(line[9]);
			CT = Double.valueOf(line[10]);
			GG = Double.valueOf(line[11]);
			GT = Double.valueOf(line[12]);
			TT = Double.valueOf(line[13]);

			/**Count Alleles:**/
			A = (2.0*AA)+AC+AG+AT;
			C = (2.0*CC)+AC+CG+CT;
			G = (2.0*GG)+AG+CG+GT;
			T = (2.0*TT)+AT+CT+GT;
			
			/**Sum of genotype counts == sum of individuals**/
			sum = A+C+G+T;
			
			if(A == sum | C == sum | G == sum | T == sum) continue;
			
			/**Initialization for the allele frequencies:**/
			ArrayList<Double> init = new ArrayList<Double>();
			for(int i = 0;i<=3;i++) init.add(0.0);
			alleles.put(id,  init);
			
			if(A != 0) {
				alleles.get(id).set(0,A/sum);
				if(A == sum) continue;
			}
			
			if(C != 0) {
				alleles.get(id).set(1,C/sum);
				
				if((A+C) == sum) continue;
			}
			
			if(G != 0) {
				alleles.get(id).set(2,G/sum);
				if((A+C+G) == sum) continue;
			}
			
			if(T != 0) {
				alleles.get(id).set(3,T/sum);
			}
		}
		in.close();
		return alleles;
	}
	
	public void read_vcf(	String[] list,
							String dir,
							HashMap<String,ArrayList<String>> info,
							HashMap<String, ArrayList<String>> genotypes,
							HashMap<String, ArrayList<Double>> qualities,
							HashMap<String, String> reference_genotype,
							HashMap<Integer,Integer> hetero_on_chrX,
							HashMap<Integer,Integer> homo_on_chrX,
							double add) throws Exception{
		
		for(String vcf:list){
			System.out.print("read "+vcf+"\n");
			
			/**-------------------------------------**/
			/**Initialize attributes for vcf file**/
			
			HashMap<String, ArrayList<String>> current = new HashMap<String, ArrayList<String>>();
			ArrayList<String> names = new ArrayList<String>();
			string_split = "-";
			
			String this_genotype[] = null;
			String this_chromsome = "";
			String this_position = "";

			String this_reference_allele = "";
			String key = "";
			int number_of_sample = info.get(NAMES).size();

			/**-------------------------------------**/
			/**Parse vcf file**/
			
			FileReader file = null;
			try {
				file = new FileReader(dir+vcf);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			BufferedReader in = new BufferedReader(file);
			
			String row = null;
			while((row = in.readLine()) != null){
				
				if(row.startsWith("##"))	continue;
				if(row.startsWith("#CHROM")){
					String[] first_line = row.split("\t");
					for(int i = 9;i<first_line.length;i++) names.add(first_line[i]);
				}else{				

					/**Parse Current VCF Line**/
					current = get_line(row);

					if(current.isEmpty()) continue;
					
					this_chromsome = current.get(CHR).get(0);
					this_position = current.get(POS).get(0);
					this_reference_allele = current.get(REF).get(0);
					
					/**Skip sex chromosomes and mitochondrial chromosomes:**/
					if(this_chromsome.equals("Y") || this_chromsome.equals("M")) continue;
						
					/**Skip INDELS:**/
					if(skip_indels(current.get(REF)) == true) continue;
					if(skip_indels(current.get(ALT)) == true) continue;
					
					/**Check Gender: (count heterozygous calls on chrX:)**/
					if(this_chromsome.equals("X")){
						for(int i = (number_of_sample);i<(number_of_sample+current.get(FORMAT_GT).size());i++){
							
								/**
								 * Get current index:
								 */
								int j = i-number_of_sample;
								
								if(!homo_on_chrX.containsKey(i)) homo_on_chrX.put(i,0);
								if(!hetero_on_chrX.containsKey(i)) hetero_on_chrX.put(i,0);
																
								/**
								 * PAR1 and PAR 2:
								 * chrX:10,000-2,781,479 and chrX:155,701,382-156,030,895**/
								if(Integer.valueOf(this_position) >= 10001 && Integer.valueOf(this_position)<= 2781479) continue;
								if(Integer.valueOf(this_position) >= 155701383 && Integer.valueOf(this_position)<= 156030895) continue;

								this_genotype = current.get(GENOTYPE).get(j).split(",");
								
								if(this_genotype[0].equals(this_genotype[1])){
									if(this_genotype[0].equals(".*\\..*")) 				continue;	//don't count .|.
									if(this_genotype[0].equals("nc")) 					continue;	//don't count .|.
									if(this_genotype[0].equals(this_reference_allele)) 	continue;	//don't count 0|0
									homo_on_chrX.put(i, (homo_on_chrX.get(i)+1));
								}
								else hetero_on_chrX.put(i, (hetero_on_chrX.get(i)+1));
						}
					}
					
					/**Add Genotype Information to Hash:**/
					key = this_chromsome+"\t"+this_position;
					//System.out.print(key);
					if(!genotypes.containsKey(key)) genotypes.put(key, new ArrayList<String>());
					if(!qualities.containsKey(key)) qualities.put(key, new ArrayList<Double>());
					
					/**Add reference genotype (== 0|0):**/
					if(!reference_genotype.containsKey(key)) reference_genotype.put(key, this_reference_allele+","+this_reference_allele);
					
					/**Fill up Genotypes for individuals which were not present in earlier vcf file:**/
					int fill_up_number = number_of_sample-genotypes.get(key).size();
					
					if(fill_up_number>0){
						for(int j=0; j<fill_up_number;j++)	{
							genotypes.get(key).add(".,.");
							qualities.get(key).add(add);
//							genotypes.get(key).add(this_reference_allele+","+this_reference_allele);
						}
					}

					/**Fill up genotypes:**/
					for(int i = 0;i<current.get(FORMAT_GT).size();i++){
						genotypes.get(key).add(current.get(GENOTYPE).get(i));
						
						if(current.get(FORMAT_GQ).get(i).equals("NA")|current.get(FORMAT_GQ).get(i).equals(".")){
							qualities.get(key).add(add);
						}else{
							qualities.get(key).add(Math.pow(10, (-1*Double.valueOf(current.get(FORMAT_GQ).get(i))/10)));
						}
					}
				}
			}
				
			/**-------------------------------------**/
			/**Collect individual names**/
			if(names.isEmpty()){
				for(int i = 0;i<current.get(FORMAT_GT).size();i++) info.get(NAMES).add("Sample"+(info.get(NAMES).size()+1));
			}else info.get(NAMES).addAll(names);
		}
}
	
public HashMap<String, ArrayList<String>> get_line(String row) throws Exception {
		
		HashMap<String, ArrayList<String>> line = new HashMap<String, ArrayList<String>>();
		if(row.startsWith("#"))	return line;	
		String[] current_row =  row.split("\t");
		if(current_row.length<9)  throw new Exception("VCF file has unexpected number of lines.\n");
			
		String[] format = null;
		
		/**initialize:**/
		line.put(CHR, new ArrayList<String>());
		line.put(POS, new ArrayList<String>());
		line.put(REF, new ArrayList<String>());
		line.put(ALT, new ArrayList<String>());
		line.put(GENOTYPE, new ArrayList<String>());
		
		/**-------------------------------------**/
		/**Define fixed Attributes of a VCF file**/
	
		line.get(CHR).add(current_row[0].replaceAll("chr|Chr", ""));
		line.get(POS).add(current_row[1]);
		for(String e: current_row[3].split(",")) line.get(REF).add(e.toUpperCase()); 
		for(String e: current_row[4].split(",")) line.get(ALT).add(e.toUpperCase()); 
				
		/**-------------------------------------**/
		/**Define variable Attributes of a VCF file**/
		
		/**FORMAT column:**/
		if(current_row[8].contains(":")){
			format = current_row[8].split(":");
			
			for(String e: format)	line.put(("FORMAT_"+e), new ArrayList<String>());
			
			/**Individual column:**/
			for(int i = 9;i<current_row.length;i++){
				String[] ind = current_row[i].split(":");
				
				/**If some individuals do not contain more information:**/
				if(ind.length < format.length){
					for(String a:format) {
						if(a.equals("GT")) continue;
						line.get("FORMAT_"+a).add("NA");
					}
				}
				
				for(int j = 0;j<ind.length;j++) {
					for(String a:ind[j].split(":")) {
						line.get("FORMAT_"+format[j]).add(a);
					}
				}
			}
			
			/** add missing formats:**/
			if(!format.toString().contains("GQ")){
				line.put(("FORMAT_"+"GQ"), new ArrayList<String>());
				
				for(int i = 9;i<current_row.length;i++){
					line.get("FORMAT_"+"GQ").add("NA");
				}
			}
			if(!format.toString().contains("DP")){
				line.put(("FORMAT_"+"DP"), new ArrayList<String>());
				
				for(int i = 9;i<current_row.length;i++){
					line.get("FORMAT_"+"DP").add("NA");
				}
			}
			if(!format.toString().contains("AD")){
				line.put(("FORMAT_"+"AD"), new ArrayList<String>());
				
				for(int i = 9;i<current_row.length;i++){
					line.get("FORMAT_"+"AD").add("NA");
				}
			}
			
		}else {
			line.put(("FORMAT_GT"), new ArrayList<String>());
			
			/**Individual column:**/
			for(int i = 9;i<current_row.length;i++){
				String[] ind = current_row[i].split(":");
				
				line.get("FORMAT_GT").add(ind[0]);
			}
		}
		
		/**Get Actual Genotype:**/
		int count = 0;
		for(String a:line.get(FORMAT_GT)){
			
			/**Check for Depth:
			 * Set to ./. == not covered, if coverage below 'minimum_depth'.
			 */

			/**Use DP flag = overall per base coverage:**/
			if(line.containsKey(FORMAT_DP) && !line.get(FORMAT_DP).get(count).matches("NA")){
				if(line.get(FORMAT_DP).get(count).equals(".") | line.get(FORMAT_DP).get(count).equals("NA")) a = ".|.";
				else if(Integer.valueOf(line.get(FORMAT_DP).get(count)) <= minimum_depth){
					a = ".|.";				
				}
			}
			
			/**Use AD flag = overall per base coverage:**/
			if(line.containsKey(FORMAT_AD) && !line.get(FORMAT_AD).get(count).matches("NA")){
				if(line.get(FORMAT_AD).get(count).equals(".") | line.get(FORMAT_AD).get(count).equals("NA")) a = ".|.";
				else {
					String[] FORMAT_AD_single = line.get(FORMAT_AD).get(count).split(",");
					if(!a.matches(".*\\..*")){
						String[] this_genotype = a.split("\\||\\/");
						if(this_genotype.length==1) {
							a = ".|.";	
							this_genotype = a.split("\\||\\/");
						}
						if(!this_genotype[0].equals(this_genotype[1])){
							if(Integer.valueOf(FORMAT_AD_single[0]) < minimum_depth | Integer.valueOf(FORMAT_AD_single[1]) < minimum_depth){
								a = ".|.";				
							}
						}
					}
				}
			}
			
			count++;
			
			/**Get string split pattern: **/
			if(a.matches(".*\\|.*")) string_split = "\\|";
			else if (a.matches(".*\\/.*"))  string_split = "\\/";
			else{
				
				/**In hemizygote chrX stat:**/
				a = a+"|"+a;
				string_split = "\\|";
			}

			List<String> tmp = new ArrayList<String>();
			for(String e:a.split(string_split)){
				if(e.equals("0")) tmp.add(line.get(REF).get(0));
				else if(e.contains(".")) tmp.add("nc");						//not covered
//				else if(e.contains(".")) tmp.add(".");						//not covered
				else tmp.add(line.get(ALT).get(Integer.valueOf(e)-1));
			}

			Collections.sort(tmp);
			line.get(GENOTYPE).add(tmp.get(0)+","+tmp.get(1));			
		}
		return line;
	}

	void update_genotypes(	HashMap<String, ArrayList<String>> genotypes,
							HashMap<String,ArrayList<String>> info,
							HashMap<String, String> reference_genotype,
							HashMap<String, ArrayList<Double>> qualities,
							double add){
	
		int number_of_individuals = info.get("NAMES").size();
			for(String key:genotypes.keySet()){
				if(genotypes.get(key).size() != number_of_individuals){
					for(int i = genotypes.get(key).size(); i<number_of_individuals;i++){
//						genotypes.get(key).add(reference_genotype.get(id));
						genotypes.get(key).add(".,.");	
						qualities.get(key).add(add);	
				}
			}
		}
	}
	
	Boolean skip_indels(ArrayList<String> str){
		Boolean skip = false;
		
		for(String e:str) if(	e.matches(".*[A-Za-z]{2,}.*") |
								e.equals("-") |
								e.equals("N") |
								e.equals(".") |
								e.equals("*")
							) skip = true;
				
		return(skip);
	}
}
   