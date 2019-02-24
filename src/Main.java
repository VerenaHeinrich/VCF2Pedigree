
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.HashMap;

public class Main{
	public static String in;												
	public static String out;												

	/**-----------------------------------------------**/
	/**Initialize additional information for vcf file:**/
	private static String NAMES = "NAMES";
	private static String HET_RATIO_CHRX = "HET_RATIO_CHRX";
	
	public static void main(String[] args) throws Exception{
		
		/**-***************************************-**/
		/**Start Time Count:**/
			long startTime = System.nanoTime();
			System.out.println("\n------------------------------------------------------------------------------------");
			
		/**-***************************************-**/
		/**Constructors:
		/**-***************************************-**/
			Commands com = new Commands();	
			Read r = new Read();
			LogLikelihood lr = new LogLikelihood();
			PedigreeBuilder pb = new PedigreeBuilder();
			Out o = new Out();
			
		/**-***************************************-**/
		/**Command Line Options:
		/**-***************************************-**/
			ArrayList<String> commands = com.get_commands(args);
			String dir  = commands.get(0);				//directory to vcf files
			final String out = commands.get(1);			//output name
			final String gt_count = commands.get(2);	//file with genotype counts for each position in the 1KG data
			
			/**Get list of vcf files from specified directory:**/
			String[] vcf_list = "".split("");
			
			if(dir.endsWith(".vcf")){
				String[] tmp = dir.split("/");
				
				dir = "";
				for(String e:tmp){
					if(!e.endsWith(".vcf")) dir = dir+e+"/";
					else{
						String[] vcf_list_tmp = e.split("/");
						vcf_list = vcf_list_tmp;
					}
				}				
			}else{
			
				/**Get list of vcf files from specified directory:**/
				vcf_list = get_file_names(dir);
			}
		
		/**Certain Error Rate for 'wrong' genotype combinations in a model:**/
			double add = 0.001;
			
		/**-*************************************************-**/
		/**Initialize LR output matrices for all three models:
		/**-*************************************************-**/

			ArrayList<ArrayList<Double>> LR_halfsib = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> LR_parent_child = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> LR_sibling = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> LR_replicate = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> LR_parent_child_sibling = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> LR_replicates_siblings = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> LR_fullSib_halfSib = new ArrayList<ArrayList<Double>>();
		
			/**Initialize normalizing matrices:**/
			ArrayList<ArrayList<Double>> Count_halfsib = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> Count_parent_child = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> Count_sibling = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> Count_replicate = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> Count_parent_child_sibling = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> Count_replicates_siblings = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> Count_fullSib_halfSib = new ArrayList<ArrayList<Double>>();
			
		/**-------------------------------------------------------------------------------**/
		/**Initialize Hash with Reference Genotypes:**/
			HashMap<String, String> reference_genotype = new HashMap<String, String>(); 
			
		/**-***************************************-**/
		/**Read genotype counts + get allele frequencies:
		/**-***************************************-**/
			 HashMap<String, ArrayList<Double>> alleles = r.read_GTCounts(	gt_count,
					 														reference_genotype); 

		/**-***************************************-**/
		/**Attributes of INFO-hash:**/
			HashMap<String,ArrayList<String>> info = new HashMap<String, ArrayList<String>>();
			info.put(NAMES, new ArrayList<String>());
			info.put(HET_RATIO_CHRX, new ArrayList<String>());

			
		/**-------------------------------------------------------------------------------**/
		/**Initialize Hash with heterozygous calls on chrX (Gender Estimation:):**/
			HashMap<Integer,Integer> hetero_on_chrX = new HashMap<Integer, Integer>();
			HashMap<Integer,Integer> homo_on_chrX = new HashMap<Integer, Integer>();

		/**-***************************************-**/
		/**Read genotype counts + get allele frequencies:
		/**-***************************************-**/
			HashMap<String, ArrayList<String>> genotypes = new HashMap<String, ArrayList<String>>();
			HashMap<String, ArrayList<Double>> qualities = new HashMap<String, ArrayList<Double>>();
			
			r.read_vcf(	vcf_list, 
						dir, 
						info,
						genotypes, 
						qualities,
						reference_genotype, 
						hetero_on_chrX, 
						homo_on_chrX,
						add);
			
		/**-***************************************-**/
		/**Test for number of genotype entries:
		/**-***************************************-**/
			if(genotypes.size()<=10000){
				System.out.print("Warning: vcf file has less than 10000 rows ("+genotypes.size()+").\n\n");
			}

		/**-------------------------------------------------------------------------------**/
		/**Update genotypes: (positions which are not covered in some individuals)**/
			r.update_genotypes(	genotypes, 
								info, 
								reference_genotype,
								qualities,
								add);
		
			
		/**-***************************************-**/
		/**LogLikelihood genotypes and calculate LOD scores:
		/**-***************************************-**/		
			
			lr.get_loglikelihood(	genotypes,
									reference_genotype,
									alleles, 
									info.get(NAMES).size(),
									LR_halfsib,
									LR_parent_child,
									LR_sibling,
									LR_replicate,
									LR_parent_child_sibling,
									LR_replicates_siblings,
									LR_fullSib_halfSib,
									Count_halfsib,
									Count_parent_child,
									Count_sibling,
									Count_replicate,
									Count_parent_child_sibling,
									Count_replicates_siblings,
									Count_fullSib_halfSib,
									add,
									qualities);
			
			/**
			 * Output:
			 */
			System.out.print("Output: HalfSibling Model"+"\n");
			o.print_out(LR_halfsib, info,out+"_HalfSib.txt");
			o.print_out(Count_halfsib, info,out+"_CountHalfSib.txt");
			
			System.out.print("Output: Parent Child Model"+"\n");
			o.print_out(LR_parent_child, info,out+"_ParCh.txt");
			o.print_out(Count_parent_child, info,out+"_CountParCh.txt");

			System.out.print("Output: Sibling Model"+"\n");
			o.print_out(LR_sibling, info, out+"_Sib.txt");
			o.print_out(Count_sibling, info,out+"_CountSib.txt");

			System.out.print("Output: Replicate Model"+"\n");
			o.print_out(LR_replicate, info, out+"_Rep.txt");
			o.print_out(Count_replicate, info,out+"_CountRep.txt");

			System.out.print("Output: Parent-Child vs Sibling Model"+"\n");
			o.print_out(LR_parent_child_sibling, info, out+"_ParCh_Sib.txt");
			
			System.out.print("Output: Replicate Model vs Sibling Model"+"\n");
			o.print_out(LR_replicates_siblings, info, out+"_Rep_Sib.txt");	
			
			System.out.print("Output: Full-Sibling Model ~ Half-Sibling Model"+"\n");
			o.print_out(LR_fullSib_halfSib, info, out+"_FullSib_HalfSib.txt");	
			
			/**Create Pedigree:**/
			HashMap<String, HashMap<String,ArrayList<String>>> pedigree = new HashMap<String, HashMap<String,ArrayList<String>>>();
			pb.get_relations(pedigree,
							info.get(NAMES).size(), 
							info.get(NAMES),
							LR_halfsib,
							LR_parent_child,
							LR_sibling,
							LR_replicate,
							LR_parent_child_sibling,
							LR_replicates_siblings,
							LR_fullSib_halfSib,
							hetero_on_chrX, 
							homo_on_chrX);
				
			o.print_ped(pedigree, out+"_ped.txt");

		/**-***************************************-**/
		/**Print Out Processing Time:**/
			long endTime = System.nanoTime();
			System.out.println("------------------------------------------------------------------------------------\nFinished!\nRuntime: "+(endTime - startTime)/1000000000 + " seconds"); 
	}
	
	/**-***************************************************************************************-**/
	/**Funktions:**/
	
	/**-----------------------------**/
	/**Get all vcf files in a directory
	 * @throws Exception **/
	public static String[] get_file_names(String dir) throws Exception{
			
		File file = new File(dir);
		if(!file.isDirectory()) throw new Exception("Directory "+dir+"does not exist.\n");
		
		String[] list = file.list(new FilenameFilter() {
			public boolean accept(File d, String name) {	return (name.endsWith(".vcf")|name.endsWith(".vcf.txt"));	}
		});	
		
		java.util.Arrays.sort(list);
		
		return list;
	}
}