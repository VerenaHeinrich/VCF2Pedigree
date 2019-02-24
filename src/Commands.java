import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


@SuppressWarnings("deprecation")
public class Commands{
	public static String in = "";
	public static String out = "";
	public static String gt_count = "";
	
	public ArrayList<String> get_commands(String[] args) throws IOException{
		try {
		      Options opt = new Options();
		      opt.addOption("d", "dir", true, "directory to vcf files or a single vcf file");
		      opt.addOption("o", "out", true, "output file name");
		      opt.addOption("g", "gt count", true, "file with gt counts from 1KG data");

		      GnuParser parser = new GnuParser();
		      CommandLine cl = parser.parse(opt, args);
		      
		      if(cl.hasOption('d')){in = cl.getOptionValue('d');}
		      else{
			         HelpFormatter f = new HelpFormatter();
			         f.printHelp("VCF2Pedigree.jar", opt);
			         System.exit(0);
			  }
		     
		      if(cl.hasOption('g')){gt_count = cl.getOptionValue('g');}
		      else{
			         HelpFormatter f = new HelpFormatter();
			         f.printHelp("VCF2Pedigree.jar", opt);
			         System.exit(0);
			  }
		      
		      if(cl.hasOption('o')){out = cl.getOptionValue('o');}

		      
		}catch (ParseException e) {
		      e.printStackTrace();
		}
		
		ArrayList<String> ret = new ArrayList<String>();
		ret.add(in);
		ret.add(out);
		ret.add(gt_count);

		return ret;
	}
}