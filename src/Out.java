import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Out{

	/**-*************************************************_**
	 * Identifier:
	 * -*************************************************-**/
	
	private static String NAMES = "NAMES";
	private static String GENDER = "GENDER";
	private static String PARENT_CHILD = "PARENT_CHILD";
	private static String SECOND_ORDER = "SECOND_ORDER";
	private static String REPLICATES = "REPLICATES";
	private static String SIBLINGS = "SIBLINGS";
	private static String FATHER = "FATHER";
	private static String MOTHER = "MOTHER";
	private static String CHILDREN = "CHILDREN";;
	private static String RATIO = "RATIO";
	
	/**Output
	 * @throws IOException **/
	public void print_out(	ArrayList<ArrayList<Double>> LR, 
							HashMap<String,ArrayList<String>> info, 
							String out)throws InterruptedException, IOException{
		
		File f = new File(out);
		FileWriter w = new FileWriter(f);
		StringBuffer o = new StringBuffer();
		
		for(int i=0;i<info.get(NAMES).size();i++)	o.append(info.get(NAMES).get(i)+"\t");
		o.append("\n");
		
		for(int i=0;i<info.get(NAMES).size();i++){
			for(int j=0;j<info.get(NAMES).size();j++){
				o.append(LR.get(i).get(j));
				if(j<info.get(NAMES).size()) o.append("\t");
			}
			o.append("\n");
		}
		
		w.write(o.toString());
		w.close();
	}
	
	public void print_ped(	HashMap<String, HashMap<String,ArrayList<String>>> pedigree,
							String out) throws InterruptedException, IOException{
		
		if(out.equals("_ped.txt")) out = "ped.txt";

		File f = new File(out);
		FileWriter w = new FileWriter(f);
		StringBuffer o = new StringBuffer();
				
		/**Header:**/
		o.append("#ID\t" +
				"MOTHER\t" +
				"FATHER\t" +
				"CHILDREN\t" +
				"SIBLINGS\t" +
				"PARENT_CHILD_RELATIONSHIP\t"+
				"REPLICATES\t" +
				"SECOND_ORDER\t" +
//				"UNRELATED\t"+
				"GENDER\t" +
				"GENDER_RATIO"+
				"\n");
		
		for(String id : pedigree.keySet()){
			o.append(id+"\t");
			
			if(!pedigree.get(id).get(MOTHER).isEmpty()) o.append(pedigree.get(id).get(MOTHER).get(0)+"\t");
			else o.append("-\t");
			
			if(!pedigree.get(id).get(FATHER).isEmpty()) o.append(pedigree.get(id).get(FATHER).get(0)+"\t");
			else o.append("-\t");
			
			if(pedigree.get(id).get(CHILDREN).isEmpty()) o.append("-\t");
			else{
				for(int i=0;i<pedigree.get(id).get(CHILDREN).size();i++){
					o.append(pedigree.get(id).get(CHILDREN).get(i));
					if(i != (pedigree.get(id).get(CHILDREN).size()-1)) o.append(",");
					else o.append("\t");
				}
			}
			
			if(pedigree.get(id).get(SIBLINGS).isEmpty()) o.append("-\t");
			else{
				for(int i=0;i<pedigree.get(id).get(SIBLINGS).size();i++){
					o.append(pedigree.get(id).get(SIBLINGS).get(i));
					if(i != (pedigree.get(id).get(SIBLINGS).size()-1)) o.append(",");
					else o.append("\t");
				}
			}
			
			if(pedigree.get(id).get(PARENT_CHILD).isEmpty()) o.append("-\t");
			else{
				for(int i=0;i<pedigree.get(id).get(PARENT_CHILD).size();i++){
					o.append(pedigree.get(id).get(PARENT_CHILD).get(i));
					if(i != (pedigree.get(id).get(PARENT_CHILD).size()-1)) o.append(",");
					else o.append("\t");
				}
			}
			
			if(pedigree.get(id).get(REPLICATES).isEmpty()) o.append("-\t");
			else{
				for(int i=0;i<pedigree.get(id).get(REPLICATES).size();i++){
					o.append(pedigree.get(id).get(REPLICATES).get(i));
					if(i != (pedigree.get(id).get(REPLICATES).size()-1)) o.append(",");
					else o.append("\t");
				}
			}
			
			if(pedigree.get(id).get(SECOND_ORDER).isEmpty()) o.append("-\t");
			else{
				for(int i=0;i<pedigree.get(id).get(SECOND_ORDER).size();i++){
					o.append(pedigree.get(id).get(SECOND_ORDER).get(i));
					if(i != (pedigree.get(id).get(SECOND_ORDER).size()-1)) o.append(",");
					else o.append("\t");
				}
			}

			o.append(pedigree.get(id).get(GENDER).get(0)+"\t");
			o.append(Math.round(Double.valueOf(pedigree.get(id).get(RATIO).get(0))*1000)/1000.0+"\n");
		}
		
		System.out.print(o+"\n");
		
		w.write(o.toString());
		w.close();
	}
}