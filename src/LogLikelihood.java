import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class LogLikelihood {
	
	/**initialize genotypes:**/
	public static String first_gt = "";
	public static String second_gt = "";
	
	/**initialize Quality Values:**/
	public static Double first_qual = 0.0;
	public static Double second_qual = 0.0;
	
	/**Initialize allele vector for each genotype:**/
	public static ArrayList<String> current_alleles_1 = new ArrayList<String>();
	public static ArrayList<String> current_alleles_2 = new ArrayList<String>();
	
	public static int current_alleles_1_first;
	public static int current_alleles_1_second;
	
	public static int current_alleles_2_first;
	public static int current_alleles_2_second;

	/**Initialize frequencies for both alleles:**/
	public static double frequency_1_first = 0.0;
	public static double frequency_1_second = 0.0;
	public static double frequency_2_first = 0.0;
	public static double frequency_2_second = 0.0;

	/**Stores the likelihood ratios for each model:**/
	public static ArrayList<Double> LR_current = new ArrayList<Double>();
	
	/**Splitted Genoytpes:**/
	HashMap<String, ArrayList<String>> splitted_genotypes = new HashMap<String,ArrayList<String>>();
	
	public void get_loglikelihood(	HashMap<String, ArrayList<String>> genotypes,
									HashMap<String, String> reference_genotype,
									HashMap<String, ArrayList<Double>> alleles,
									int number_of_individuals,
									ArrayList<ArrayList<Double>> LR_halfsib,
									ArrayList<ArrayList<Double>> LR_parent_child,
									ArrayList<ArrayList<Double>> LR_sibling,
									ArrayList<ArrayList<Double>> LR_replicate,
									ArrayList<ArrayList<Double>> LR_parent_child_sibling,
									ArrayList<ArrayList<Double>> LR_replicates_siblings,
									ArrayList<ArrayList<Double>> LR_fullSib_halfSib,
									ArrayList<ArrayList<Double>> Count_halfsib,
									ArrayList<ArrayList<Double>> Count_parent_child,
									ArrayList<ArrayList<Double>> Count_sibling,
									ArrayList<ArrayList<Double>> Count_replicate,
									ArrayList<ArrayList<Double>> Count_parent_child_sibling,
									ArrayList<ArrayList<Double>> Count_replicates_siblings,
									ArrayList<ArrayList<Double>> Count_fullSib_halfSib,
									Double add,
									HashMap<String, ArrayList<Double>> qualities) throws IOException, InterruptedException{
		
		/**Initialize LR matrices:**/
		LR_halfsib = init(LR_halfsib,number_of_individuals);
		LR_parent_child = init(LR_parent_child,number_of_individuals);
		LR_sibling = init(LR_sibling,number_of_individuals);
		LR_replicate = init(LR_replicate,number_of_individuals);
		LR_parent_child_sibling = init(LR_parent_child_sibling,number_of_individuals);
		LR_replicates_siblings = init(LR_replicates_siblings,number_of_individuals);
		LR_fullSib_halfSib = init(LR_fullSib_halfSib,number_of_individuals);

		Count_halfsib = init(Count_halfsib,number_of_individuals);
		Count_parent_child = init(Count_parent_child,number_of_individuals);
		Count_sibling = init(Count_sibling,number_of_individuals);
		Count_replicate = init(Count_replicate,number_of_individuals);
		Count_parent_child_sibling = init(Count_parent_child_sibling,number_of_individuals);
		Count_replicates_siblings = init(Count_replicates_siblings,number_of_individuals);
		Count_fullSib_halfSib = init(Count_fullSib_halfSib,number_of_individuals);
		
		/**Define allele key:**/
		HashMap<String, Integer> allele_key = new HashMap<String, Integer>();
		allele_key.put("A", 0);
		allele_key.put("C", 1);
		allele_key.put("G", 2);
		allele_key.put("T", 3);

		/**Get splittet genoytpes**/
		for(String s1: allele_key.keySet()){
			for(String s2: allele_key.keySet()){
				splitted_genotypes.put(s1+","+s2, new ArrayList<String>());
				splitted_genotypes.get(s1+","+s2).add(s1);
				splitted_genotypes.get(s1+","+s2).add(s2);
			}
		}

		/**Go through every position two compare each pair of genotypes:**/
		for(String pos:genotypes.keySet()){

			/**get allele frequencies at this position:
			 * 
			 * NOTE: Just for positions in genotype vector!
			 * 
			 * **/
			if(!alleles.containsKey(pos)) {				
				continue;
			}
	
			ArrayList<Double> current_allele_frequencies = alleles.get(pos);
		
			/**Check every pair of genoytpes:**/
			for(int i=0;i<=(number_of_individuals-1);i++){
				for(int j=(i+1);j<number_of_individuals;j++){
					
					/**Get genotypes for comparison:**/
					first_gt = genotypes.get(pos).get(i);
					second_gt = genotypes.get(pos).get(j);
						
					/**skip comparisons where at least one genotype is not covered:**/
					if(first_gt.matches(".*\\.,\\..*") & second_gt.matches(".*\\.,\\..*")) continue;
					if(first_gt.matches(".*nc,nc.*") | second_gt.matches(".*nc,nc.*")) continue;

					if(first_gt.matches(".*\\.,\\..*")){
						first_gt = reference_genotype.get(pos);	
					}
					if(second_gt.matches(".*\\.,\\..*")){
						second_gt = reference_genotype.get(pos);
					}

					/**Add Quality Value:**/
					first_qual = qualities.get(pos).get(i);
					second_qual = qualities.get(pos).get(j);
					
					/**Take mean quality score**/
					add=0.001;//(first_qual+second_qual)/2;
					
					/**allele vector for each genotype:**/
					current_alleles_1 = splitted_genotypes.get(first_gt);
					current_alleles_2 = splitted_genotypes.get(second_gt);
				
					current_alleles_1_first = allele_key.get(current_alleles_1.get(0));
					current_alleles_1_second = allele_key.get(current_alleles_1.get(1));
					
					current_alleles_2_first = allele_key.get(current_alleles_2.get(0));
					current_alleles_2_second = allele_key.get(current_alleles_2.get(1));

					frequency_1_first = current_allele_frequencies.get(current_alleles_1_first);
					frequency_1_second = current_allele_frequencies.get(current_alleles_1_second);

					frequency_2_first = current_allele_frequencies.get(current_alleles_2_first);
					frequency_2_second = current_allele_frequencies.get(current_alleles_2_second);

					/**check if homozygous/heterozygous genotypes:**/
					if(current_alleles_1_first == current_alleles_2_first & current_alleles_1_second == current_alleles_2_second){
						
						/**both genotypes are equal homozygous:**/
						if(current_alleles_1_first == current_alleles_1_second){	
							LR_current=compare_hom_hom_equal(frequency_1_first, frequency_1_first);
						}
						
						/**both genotypes are heterozygous:**/
						else{
							LR_current=compare_het_het(frequency_1_first, frequency_1_second, add);
						}
					}else{

						/**first genotype is homozygous**/
						if(current_alleles_1_first == current_alleles_1_second){
							
							/**second genotype is homozygous**/
							if(current_alleles_2_first == current_alleles_2_second){	
								LR_current=compare_hom_hom_unequal(frequency_1_first, frequency_2_first, add);
							}
							
							/**one allele of second genotype equals both alleles of first genotype:**/
							else{
								if(current_alleles_2_first == current_alleles_1_first){		
									LR_current=compare_hom_het(frequency_1_first, frequency_2_second, add);
								}else if(current_alleles_2_second == current_alleles_1_first){
									LR_current=compare_hom_het(frequency_1_first, frequency_2_first, add);
								}else{
									
									/**no allele of second genotype equals allele of first genotype:**/
									LR_current=compare_hom_het_unequal(current_alleles_1_first, current_alleles_2_first, current_alleles_2_second, add);
								}
							}
						}
						/**second genotype is homozygous**/
						else if(current_alleles_2_first == current_alleles_2_second){
							
							/**first genotype is homozygous**/
							if(current_alleles_1_first == current_alleles_1_second){								
								LR_current=compare_hom_hom_unequal(frequency_2_first, frequency_1_first, add);
							}

							/**one allele of second genotype equals both allels of first genotype:**/
							else{
								if(current_alleles_1_first == current_alleles_2_first){									
									LR_current=compare_hom_het(frequency_2_first, frequency_1_second, add);
								}else if(current_alleles_1_second == current_alleles_2_first){
									LR_current=compare_hom_het(frequency_2_first, frequency_1_first, add);
								}else{
									
									/**no allele of first genotype equals allele of second genotype:**/
									LR_current=compare_hom_het_unequal(current_alleles_2_first, current_alleles_1_first, current_alleles_1_second, add);
								}
							}
						}else{
							
							/**Both Genotypes are heterozygous but not equal:**/							
							if(current_alleles_1_first != current_alleles_2_first & current_alleles_1_first != current_alleles_2_second){
								
								/**First Case: all alleles are different:**/
								if(current_alleles_1_second != current_alleles_2_first & current_alleles_1_second != current_alleles_2_second){
									LR_current=compare_het_het_all_unequal(current_alleles_1_first, current_alleles_1_second,current_alleles_2_first, current_alleles_2_second, add);
								}else{
									
									/**one allele is equal in both genotypes:**/
									if(current_alleles_1_second == current_alleles_2_first){
										LR_current=compare_het_het_one_unequal(current_alleles_1_second, current_alleles_1_first,current_alleles_2_second, add);
									}else if(current_alleles_1_second == current_alleles_2_second){
										LR_current=compare_het_het_one_unequal(current_alleles_1_second, current_alleles_1_first,current_alleles_2_first, add);
									}
								}
							}else{
								if(current_alleles_1_first == current_alleles_2_first){
									LR_current=compare_het_het_one_unequal(current_alleles_1_first, current_alleles_1_second,current_alleles_2_second, add);
								}else if(current_alleles_1_first == current_alleles_2_second){
									LR_current=compare_het_het_one_unequal(current_alleles_1_first, current_alleles_1_second,current_alleles_2_first, add);
								}
							}
						}
					}
				
					if(LR_current.isEmpty()) continue;
					
					/**Count positions for normalizing:**/
					Count_halfsib.get(i).set(j, Count_halfsib.get(i).get(j)+1.0);
					Count_parent_child.get(i).set(j, Count_parent_child.get(i).get(j)+1.0);
					Count_sibling.get(i).set(j, Count_sibling.get(i).get(j)+1.0);
					Count_replicate.get(i).set(j, Count_replicate.get(i).get(j)+1.0);
					Count_parent_child_sibling.get(i).set(j, Count_parent_child_sibling.get(i).get(j)+1.0);
					Count_replicates_siblings.get(i).set(j, Count_replicates_siblings.get(i).get(j)+1.0);
					Count_fullSib_halfSib.get(i).set(j, Count_fullSib_halfSib.get(i).get(j)+1.0);
					
					/**add current LR to output matrices:**/
					LR_halfsib.get(i).set(j, LR_halfsib.get(i).get(j)+LR_current.get(0));
					LR_parent_child.get(i).set(j, LR_parent_child.get(i).get(j)+LR_current.get(1));
					LR_sibling.get(i).set(j, LR_sibling.get(i).get(j)+LR_current.get(2));
					LR_replicate.get(i).set(j, LR_replicate.get(i).get(j)+LR_current.get(3));
					LR_parent_child_sibling.get(i).set(j, LR_parent_child_sibling.get(i).get(j)+LR_current.get(4));
					LR_replicates_siblings.get(i).set(j, LR_replicates_siblings.get(i).get(j)+LR_current.get(5));
					LR_fullSib_halfSib.get(i).set(j, LR_fullSib_halfSib.get(i).get(j)+LR_current.get(6));	
					
				}
			}
		}

//		/**normalize by number of positions:**/
//		for(int i=0;i<=(number_of_individuals-1);i++){
//			for(int j=(i+1);j<number_of_individuals;j++){
//				
//				LR_halfsib.get(i).set(j, LR_halfsib.get(i).get(j)/(Count_halfsib.get(i).get(j)));
//				LR_parent_child.get(i).set(j, LR_parent_child.get(i).get(j)/(Count_parent_child.get(i).get(j)));
//				LR_sibling.get(i).set(j, LR_sibling.get(i).get(j)/(Count_sibling.get(i).get(j)));
//				LR_replicate.get(i).set(j, LR_replicate.get(i).get(j)/(Count_replicate.get(i).get(j)));
//				LR_parent_child_sibling.get(i).set(j, LR_parent_child_sibling.get(i).get(j)/(Count_parent_child_sibling.get(i).get(j)));
//				LR_replicates_siblings.get(i).set(j, LR_replicates_siblings.get(i).get(j)/(Count_replicates_siblings.get(i).get(j)));
//				LR_fullSib_halfSib.get(i).set(j, LR_fullSib_halfSib.get(i).get(j)/(Count_fullSib_halfSib.get(i).get(j)));
//
//			}
//		}
	}
	
	/**
	 * Initialize Arary
	 * @param array
	 * @param number_of_individuals
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<ArrayList<Double>> init(ArrayList<ArrayList<Double>> array, int number_of_individuals)throws InterruptedException{
		for(int i=0;i<number_of_individuals;i++){
			array.add(new ArrayList<Double>());
			
			for(int j=0;j<number_of_individuals;j++){
				array.get(i).add(0.0);
			}
		}
		return array;
	}

	/**
	 * Get Ratio for two un-equal homozygous Genotypes
	 * @param f1
	 * @param f2
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<Double> compare_hom_hom_unequal(double f1,double f2, double add)throws InterruptedException{
		ArrayList<Double> LR = new ArrayList<Double>();
		
		if(f1 == 0| f2 == 0){
			return(LR);
		}

			/**Model: Halfsib ~ Unrelated:**/
				LR.add(Math.log10(0.5));
			
			/**Model: Parent ~ Child:**/
				LR.add(Math.log10(add/(2*f1*f1*f2*f2)));
					
			/**Model: siblings:**/
				LR.add(Math.log10(0.25));
	
			/**Model: Replicates:**/
				LR.add(Math.log10(add/(f1*f1*f2*f2)));
				
			/**Model: Parent-Child ~ Siblings:**/
				LR.add(Math.log10((2*add)/(f1*f1*f2*f2)));
				
			/**Model: Replicates ~ Siblings:**/
				LR.add(Math.log10((4*add)/(f1*f1*f2*f2)));
				
			/**Model: Full Siblings ~ Half Siblings:**/
				LR.add(Math.log10(0.5));
	
			return LR;
	}
	
	/**
	 *Get Ratio for two equal homozygous Genotypes
	 * @param f1
	 * @param f2
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<Double> compare_hom_hom_equal(double f1, double f2)throws InterruptedException{
		ArrayList<Double> LR = new ArrayList<Double>();
		
		if(f1==0){
			return(LR);
		}

			/**Model: Halfsib ~ Unrelated:**/
				LR.add(Math.log10((f1+1)/(2*f1)));

			/**Model: Parent ~ Child:**/
				LR.add(Math.log10((1/f1)));
	
			/**Model: siblings:**/
				LR.add(Math.log10(((f1+1)*(f1+1)/(4*f1*f1))));
	
			/**Model: Replicates:**/
				LR.add(Math.log10((1/(f1*f1))));
				
			/**Model: Parent-Child ~ Siblings:**/
				LR.add(Math.log10((4*f1)/((f1+1)*(f1+1))));
			
			/**Model: Replicates ~ Siblings:**/
				LR.add(Math.log10(4/((f1+1)*(f1+1))));

			/**Model: Full Siblings ~ Half Siblings:**/
				LR.add(Math.log10((f1+1)/(2*f1)));

			return LR;
	}
	
	/**
	 * Get Ratio for one homozygous and one heterozygous Genotype
	 * @param f1
	 * @param f2
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<Double> compare_hom_het(double f1,double f2, double add)throws InterruptedException{
		ArrayList<Double> LR = new ArrayList<Double>();
		
		if(f1 == 0| f2 == 0){		
			return(LR);
		}
		
			/**Model: Halfsib ~ Unrelated:**/
				LR.add(Math.log10((2*f1+1)/(4*f1)));

			/**Model: Parent ~ Child:**/
				LR.add(Math.log10((1/(2*f1))));
					
			/**Model: siblings:**/
				LR.add(Math.log10(((f1+1)/(4*f1))));
		
			/**Model: Replicates:**/
				LR.add(Math.log10(add/2*f1*f1*f1*f2));
			
			/**Model: Parent-Child ~ Siblings:**/
				LR.add(Math.log10(2/(f1+1)));
			
			/**Model: Replicates ~ Siblings:**/
				LR.add(Math.log10(2*add/((f1*f1*f2)*(f1+1))));
				
			/**Model: Full Siblings ~ Half Siblings:**/
				LR.add(Math.log10((f1+1)/(2*f1+1)));

			return LR;			
	}
	
	/**
	 * Get Ratio for two heterozygous Genotypes
	 * @param f1
	 * @param f2
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<Double> compare_het_het(double f1,double f2, double add)throws InterruptedException{
		ArrayList<Double> LR = new ArrayList<Double>();
		
		if(f1 == 0| f2 == 0){
			return(LR);
		}

			/**Model: Halfsib ~ Unrelated:**/
				LR.add(Math.log10((4*f1*f2+f1+f2)/(8*f1*f2)));
			
			/**Model: Parent ~ Child:**/
				LR.add(Math.log10((1/(4*f1*f2))));
	
			/**Model: siblings:**/
				LR.add(Math.log10(((2*f1*f2+f1+f2+1)/(8*f1*f2))));
				
			/**Model: Replicates:**/
				LR.add(Math.log10((1/(2*f1*f2))));
	
			/**Model: Parent-Child ~ Siblings:**/
				LR.add(Math.log10((2/(2*f1*f2+f1+f2+1))));
				
			/**Model: Replicates ~ Siblings:**/
				LR.add(Math.log10(2/(1+f1*f2)));
				
			/**Model: Full Siblings ~ Half Siblings:**/
				LR.add(Math.log10((2*f1*f2+f1+f2+1)/(4*f1*f2+f1+f2)));

			return LR;
	}
	
	/**
	 * Get Ratio for two heterozygous Genotypes
	 * @param f1
	 * @param f2
	 * @param f3
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<Double> compare_hom_het_unequal(double f1,double f2, double f3, double add)throws InterruptedException{
		ArrayList<Double> LR = new ArrayList<Double>();
		
		if(f1 == 0| f2 == 0 | f3 == 0){
			return(LR);
		}

			/**Model: Halfsib ~ Unrelated:**/
				LR.add(Math.log10(0.5));
			
			/**Model: Parent ~ Child:**/
				LR.add(Math.log10((add/(4*f1*f1*f2*f3))));
	
			/**Model: siblings:**/
				LR.add(Math.log10((0.25)));
				
			/**Model: Replicates:**/
				LR.add(Math.log10((add/(4*f1*f1*f2*f3))));
	
			/**Model: Parent-Child ~ Siblings:**/
				LR.add(Math.log10((add/(f1*f1*f2*f3))));
				
			/**Model: Replicates ~ Siblings:**/
				LR.add(Math.log10((add/(f1*f1*f2*f3))));
				
			/**Model: Full Siblings ~ Half Siblings:**/
				LR.add(Math.log10((0.5)));

			return LR;
	}
	
	/**
	 * Get Ratio for two heterozygous Genotypes
	 * @param f1
	 * @param f2
	 * @param f3
	 * @param f4
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<Double> compare_het_het_all_unequal(double f1,double f2, double f3, double f4, double add)throws InterruptedException{
		ArrayList<Double> LR = new ArrayList<Double>();
		
		if(f1 == 0| f2 == 0 | f3 == 0 | f4 == 0){
			return(LR);
		}

			/**Model: Halfsib ~ Unrelated:**/
				LR.add(Math.log10(0.5));
			
			/**Model: Parent ~ Child:**/
				LR.add(Math.log10((add/(8*f1*f2*f3*f4))));
	
			/**Model: siblings:**/
				LR.add(Math.log10((0.25)));
				
			/**Model: Replicates:**/
				LR.add(Math.log10((add/(8*f1*f2*f3*f4))));
	
			/**Model: Parent-Child ~ Siblings:**/
				LR.add(Math.log10((add/(2*f1*f2*f3*f4))));
				
			/**Model: Replicates ~ Siblings:**/
				LR.add(Math.log10((add/(2*f1*f2*f3*f4))));
				
			/**Model: Full Siblings ~ Half Siblings:**/
				LR.add(Math.log10((0.5)));

			return LR;
	}
	
	/**
	 * Get Ratio for two heterozygous Genotypes
	 * @param f1
	 * @param f2
	 * @param f3
	 * @return
	 * @throws InterruptedException
	 */
	public ArrayList<Double> compare_het_het_one_unequal(double f1,double f2, double f3, double add)throws InterruptedException{
		ArrayList<Double> LR = new ArrayList<Double>();
		
		if(f1 == 0| f2 == 0 | f3 == 0){	
			return(LR);
		}

			/**Model: Halfsib ~ Unrelated:**/
				LR.add(Math.log10((4*f1+1)/(8*f1)));
			
			/**Model: Parent ~ Child:**/
				LR.add(Math.log10((1/(8*f1))));
	
			/**Model: siblings:**/
				LR.add(Math.log10((2*f1+1)/(8*f1)));
				
			/**Model: Replicates:**/
				LR.add(Math.log10((add/(8*f1*f1*f2*f3))));
	
			/**Model: Parent-Child ~ Siblings:**/
				LR.add(Math.log10((1/(2*f1+1))));
				
			/**Model: Replicates ~ Siblings:**/
				LR.add(Math.log10((add/(f1*f2*f3*(2*f1+1)))));
				
			/**Model: Full Siblings ~ Half Siblings:**/
				LR.add(Math.log10((2*f1+1)/(4*f1+1)));

			return LR;
	}
}