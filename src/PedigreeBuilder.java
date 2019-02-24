import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class PedigreeBuilder {
	
	/**-*************************************************_**
	 * Identifier:
	 * -*************************************************-**/
	private static String GENDER = "GENDER";
	private static String PARENT_CHILD = "PARENT_CHILD";
	private static String SECOND_ORDER = "SECOND_ORDER";
	private static String REPLICATES = "REPLICATES";
	private static String SIBLINGS = "SIBLINGS";
	private static String UNRELATED = "UNRELATED";
	private static String FATHER = "FATHER";
	private static String MOTHER = "MOTHER";
	private static String CHILDREN = "CHILDREN";
	private static String PARENT_UNKNOWN_GENDER = "PARENT_UNKNOWN_GENDER";
	private static String AUNT = "AUNT";
	private static String UNCLE = "UNCLE";
	private static String MALE = "male";
	private static String FEMALE = "female";
	private static String UNKNOWN = "unknown";
	private static String RATIO = "RATIO";
	
	private String id = "";
	private String sec_id = "";
	
	private HashMap<String, String> gender_key = new HashMap<String, String>();

	/**
	 * 
	 * @param pedigree
	 * @param number_of_individuals
	 * @param names
	 * @param LR_parent_child
	 * @param LR_sibling
	 * @param LR_replicate
	 * @param LR_parent_child_sibling
	 * @param LR_replicates_siblings
	 * @param LR_fullSib_halfSib
	 * @param hetero_on_chrX
	 * @param homo_on_chrX
	 * @throws InterruptedException
	 */
	public void get_relations(	HashMap<String, HashMap<String,ArrayList<String>>> pedigree,
								int number_of_individuals,
								ArrayList<String> names,
								ArrayList<ArrayList<Double>> LR_halfsib,
								ArrayList<ArrayList<Double>> LR_parent_child,
								ArrayList<ArrayList<Double>> LR_sibling,
								ArrayList<ArrayList<Double>> LR_replicate,
								ArrayList<ArrayList<Double>> LR_parent_child_sibling,
								ArrayList<ArrayList<Double>> LR_replicates_siblings,
								ArrayList<ArrayList<Double>> LR_fullSib_halfSib,
								HashMap<Integer,Integer> hetero_on_chrX,
								HashMap<Integer,Integer> homo_on_chrX
								) throws InterruptedException{
		
		/**
		 * Define gender key:
		 */
		gender_key.put(FATHER, "male");
		gender_key.put(MOTHER, "female");
		
		/**
		 * Get Relationship degrees:
		 */
		
		for(int ind=0;ind<number_of_individuals;ind++){
			id = names.get(ind);
			
			/**Initialize Pedigree Attributes:**/
			if(!pedigree.containsKey(id)) init_pedigree_at_id(pedigree, id);
			
			/**Gender:**/
			get_gender(	hetero_on_chrX, 
						homo_on_chrX,
						pedigree,
						id,
						ind);
			
			
			/**Get Relations:**/
			for(int sec_ind=ind;sec_ind<number_of_individuals;sec_ind++){
				if(ind==sec_ind) continue;
				sec_id = names.get(sec_ind);

				/**Initialize Pedigree Attributes:**/
				if(!pedigree.containsKey(sec_id)) init_pedigree_at_id(pedigree, sec_id);
				
				/**Unrelated:**/
				if( LR_sibling.get(ind).get(sec_ind)<0 & LR_parent_child.get(ind).get(sec_ind)<0 & LR_halfsib.get(ind).get(sec_ind) <0) {
					if (! pedigree.get(id).get(UNRELATED).contains(sec_id)) pedigree.get(id).get(UNRELATED).add(sec_id);
					if (! pedigree.get(sec_id).get(UNRELATED).contains(id)) pedigree.get(sec_id).get(UNRELATED).add(id);
					continue;
				}
				
				/**Parent ~ Child**/
				if(LR_parent_child_sibling.get(ind).get(sec_ind) >= 0.0) {
					if (! pedigree.get(id).get(PARENT_CHILD).contains(sec_id)) pedigree.get(id).get(PARENT_CHILD).add(sec_id);
					if (! pedigree.get(sec_id).get(PARENT_CHILD).contains(id)) pedigree.get(sec_id).get(PARENT_CHILD).add(id);
					continue;
				}

				/**Replicates ~ Siblings:**/
				if(LR_replicates_siblings.get(ind).get(sec_ind) >= 0.0) {
					if (! pedigree.get(id).get(REPLICATES).contains(sec_id)) pedigree.get(id).get(REPLICATES).add(sec_id);
					if (! pedigree.get(sec_id).get(REPLICATES).contains(id)) pedigree.get(sec_id).get(REPLICATES).add(id);
					continue;
				}


				/**Full Siblings ~ Half Siblings:**/
				if(LR_sibling.get(ind).get(sec_ind) >= 0.0 & LR_fullSib_halfSib.get(ind).get(sec_ind) >= 0.0) {

					if (! pedigree.get(id).get(SIBLINGS).contains(sec_id)) pedigree.get(id).get(SIBLINGS).add(sec_id);
					if (! pedigree.get(sec_id).get(SIBLINGS).contains(id)) pedigree.get(sec_id).get(SIBLINGS).add(id);

				}else if (LR_fullSib_halfSib.get(ind).get(sec_ind)< 0.0) {
					if (! pedigree.get(id).get(SECOND_ORDER).contains(sec_id)) pedigree.get(id).get(SECOND_ORDER).add(sec_id);
					if (! pedigree.get(sec_id).get(SECOND_ORDER).contains(id)) pedigree.get(sec_id).get(SECOND_ORDER).add(id);
					continue;
				}
			}
		}

		/**
		 * Get Mother and Father and Children Pedigree Information:
		 */
		for(int ind=0;ind<(number_of_individuals);ind++){
			id = names.get(ind);
			
			/**
			 * Specify father:
			 */
			String father = get_parent(pedigree, id, FATHER);
			if(!father.isEmpty()) {
				if(!pedigree.get(id).get(FATHER).contains(father)) pedigree.get(id).get(FATHER).add(father);
				if(!pedigree.get(father).get(CHILDREN).contains(id)) pedigree.get(father).get(CHILDREN).add(id);
			}
			
			/**
			 * Specify mother:
			 */
			String mother = get_parent(pedigree, id, MOTHER);
			if(!mother.isEmpty()) {
				if(!pedigree.get(id).get(MOTHER).contains(mother)) pedigree.get(id).get(MOTHER).add(mother);
				if(!pedigree.get(mother).get(CHILDREN).contains(id)) pedigree.get(mother).get(CHILDREN).add(id);
			}
		}
		
		/**
		 * Get remaining children/second order(Grandparent, Grandchildren) Pedigree Information:
		 */
		for(int ind=0;ind<(number_of_individuals);ind++){
			id = names.get(ind);
		
			/**
			 * Check for Grandparents:
			 */
			ArrayList<String> parent_ids = pedigree.get(id).get(PARENT_CHILD);
			for(String parent_id:parent_ids){
				if(pedigree.get(id).get(CHILDREN).contains(parent_id)) continue;
				
				for(String gender_grandparent:gender_key.keySet()){
					if(pedigree.get(parent_id).get(gender_grandparent).isEmpty()) continue;
					
					String grand_parent = pedigree.get(parent_id).get(gender_grandparent).get(0);
					for(String gender_parent:gender_key.keySet()){
						if(!pedigree.get(parent_id).get(GENDER).get(0).equals(gender_key.get(gender_parent))) continue;
						
						if(!pedigree.get(id).get(gender_parent).contains(parent_id)) pedigree.get(id).get(gender_parent).add(parent_id);
						if(!pedigree.get(parent_id).get(CHILDREN).contains(id)) pedigree.get(parent_id).get(CHILDREN).add(id);
					}
					if(!pedigree.get(id).get(SECOND_ORDER).contains(grand_parent)) pedigree.get(id).get(SECOND_ORDER).add(grand_parent);
					if(!pedigree.get(grand_parent).get(SECOND_ORDER).contains(id)) pedigree.get(grand_parent).get(SECOND_ORDER).add(id);
				}
			}
		}
		
		/**
		 * Correct Sibling states for equal Parents (Case for inbred families):
		 */
		for(int ind=0;ind<(number_of_individuals);ind++){
			id = names.get(ind);
			ArrayList<String> children_ids = pedigree.get(id).get(CHILDREN);
			
			for(String child:children_ids){
				if(pedigree.get(child).get(MOTHER).isEmpty() | pedigree.get(child).get(FATHER).isEmpty()) continue;
				
				for(String sib:children_ids){
					if(child.equals(sib)) continue;

					if(pedigree.get(sib).get(MOTHER).isEmpty() | pedigree.get(sib).get(FATHER).isEmpty()) continue;

					/**
					 * If siblings share same mother AND same father:
					 */
						
					if(	pedigree.get(child).get(MOTHER).get(0).equals(pedigree.get(sib).get(MOTHER).get(0)) &
						pedigree.get(child).get(FATHER).get(0).equals(pedigree.get(sib).get(FATHER).get(0))){
								if(!pedigree.get(child).get(SIBLINGS).contains(sib)) pedigree.get(child).get(SIBLINGS).add(sib);
								if(!pedigree.get(sib).get(SIBLINGS).contains(child)) pedigree.get(sib).get(SIBLINGS).add(child);
					}
				}
			}
		}
	}
	
	/**
	 * Get Gender Estimate.
	 * @param hetero_on_chrX
	 * @param homo_on_chrX
	 * @param pedigree
	 * @param id
	 * @param ind
	 * @throws InterruptedException
	 */
	public void get_gender(		HashMap<Integer,Integer> hetero_on_chrX,
								HashMap<Integer,Integer> homo_on_chrX,
								HashMap<String, HashMap<String,ArrayList<String>>> pedigree,
								String id,
								int ind)throws InterruptedException{
		
		double threshold = 0.45;
		double gender_estimate = 0.0;
		
		if(hetero_on_chrX.containsKey(ind) && Double.valueOf(hetero_on_chrX.get(ind)) >0){
			gender_estimate = Double.valueOf(hetero_on_chrX.get(ind))/(Double.valueOf(hetero_on_chrX.get(ind))+Double.valueOf(homo_on_chrX.get(ind)));
		}
		pedigree.get(id).get(RATIO).add(String.valueOf(gender_estimate));
		
		if(gender_estimate >= threshold) pedigree.get(id).get(GENDER).add(FEMALE);
		else if (gender_estimate <= threshold & threshold>0) pedigree.get(id).get(GENDER).add(MALE);
		else pedigree.get(id).get(GENDER).add(UNKNOWN);
	}
	
	/**
	 * 
	 * @param pedigree
	 * @param id
	 * @param attribute
	 * @return parent (FATHER or MOTHER)
	 * @throws InterruptedException
	 */
	@SuppressWarnings("unchecked")
	public String get_parent(	HashMap<String, HashMap<String,ArrayList<String>>> pedigree, 
								String id,
								String attribute) throws InterruptedException{
		
		/**Initialize String for parent**/
		String parent = "";
		
		if(!pedigree.get(id).get(PARENT_CHILD).isEmpty()){
			ArrayList<String> parent_ids = pedigree.get(id).get(PARENT_CHILD);
			ArrayList<String> unrelated_parent_ids = get_subset(pedigree, parent_ids, UNRELATED);

			/**
			 * case:two parents:
			 */
			if(unrelated_parent_ids.size() == 2){
				for(String parent_id:unrelated_parent_ids){
					if(pedigree.get(parent_id).get(GENDER).get(0).equals(gender_key.get(attribute))){
						parent = parent_id;
					}
				}
			}
			
			/**
			 * case:one parent; at least one siblings:
			 */
			if(pedigree.get(id).get(PARENT_CHILD).size() == 1 & !pedigree.get(id).get(SIBLINGS).isEmpty()){
				ArrayList<String> siblings = pedigree.get(id).get(SIBLINGS);
				String parent_id = parent_ids.get(0);
				
				for(String sib:siblings) {
					if(pedigree.get(sib).get(PARENT_CHILD).contains(parent_id)){
						if(pedigree.get(parent_id).get(GENDER).get(0).equals(gender_key.get(attribute))){
							parent = parent_id;
						}
					}
				}
			}

			/**
			 * case:one parent; at least two children:
			 */
			if(pedigree.get(id).get(PARENT_CHILD).size() > 2  & get_subset(pedigree, parent_ids, SIBLINGS).size() >= 2){
				ArrayList<String> children = get_subset(pedigree, pedigree.get(id).get(PARENT_CHILD), SIBLINGS);

				if(children.size()>=2 & children.size() != pedigree.get(id).get(PARENT_CHILD).size()){
					@SuppressWarnings("rawtypes")
					Collection sublist_children = new ArrayList();
					sublist_children.addAll(children);
					
					@SuppressWarnings("rawtypes")
					Collection sublist_parent = new ArrayList();
					sublist_parent.addAll(pedigree.get(id).get(PARENT_CHILD));

					sublist_parent.removeAll(sublist_children);
					String parent_id = sublist_parent.iterator().next().toString();
					if(pedigree.get(parent_id).get(GENDER).get(0).equals(gender_key.get(attribute))){
						parent = parent_id;
					}
				}
			}
		}
		return parent;
	}
	
	/**
	 * Get unrelated ids from a set of given ids.
	 * @param pedigree
	 * @param ids
	 * @return list of unrelated individuals
	 * @throws InterruptedException
	 */
	public ArrayList<String> get_subset(	HashMap<String, HashMap<String,ArrayList<String>>> pedigree,
											ArrayList<String> ids,
											String Attribute) throws InterruptedException{
		
		ArrayList<String> subset_ids = new ArrayList<String>();
		for(String id1: ids){
			for(String id2: ids){
				if(id1.equals(id2)) continue;
				if(pedigree.get(id1).get(Attribute).contains(id2)) {
					if(!subset_ids.contains(id1)) subset_ids.add(id1);
					if(!subset_ids.contains(id2)) subset_ids.add(id2);
				}
			}
		}
		return subset_ids;
	}
	
	/**
	 * Initialize pedigree attributes.
	 * @param pedigree
	 * @param id
	 * @throws InterruptedException
	 */
	public void init_pedigree_at_id(HashMap<String, HashMap<String,ArrayList<String>>> pedigree,
									String id)throws InterruptedException{
		
		pedigree.put(id, new HashMap<String, ArrayList<String>>());
		pedigree.get(id).put(GENDER, new ArrayList<String>());
		pedigree.get(id).put(RATIO, new ArrayList<String>());
		pedigree.get(id).put(PARENT_CHILD, new ArrayList<String>());
		pedigree.get(id).put(SECOND_ORDER, new ArrayList<String>());
		pedigree.get(id).put(REPLICATES, new ArrayList<String>());
		pedigree.get(id).put(SIBLINGS, new ArrayList<String>());
		pedigree.get(id).put(FATHER, new ArrayList<String>());
		pedigree.get(id).put(MOTHER, new ArrayList<String>());
		pedigree.get(id).put(UNRELATED, new ArrayList<String>());
		pedigree.get(id).put(PARENT_UNKNOWN_GENDER, new ArrayList<String>());
		pedigree.get(id).put(CHILDREN, new ArrayList<String>());
		pedigree.get(id).put(AUNT, new ArrayList<String>());
		pedigree.get(id).put(UNCLE,new ArrayList<String>());
	}
}