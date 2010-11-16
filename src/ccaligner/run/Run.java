/*
 * $Id: Run,v 1.3 2005/04/03 19:38:21 ahmed Exp $
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package ccaligner.run;

import java.io.*;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.cli.*;

import org.biojava.bio.*;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.symbol.*;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.db.HashRichSequenceDB;
import org.biojavax.bio.seq.*;

import ccaligner.*;
import ccaligner.Alignment;
import ccaligner.matrix.*;
import ccaligner.formats.Pair;

/**
 * Example of using JAligner API to align P53 human against
 * P53 mouse using Smith-Waterman-Gotoh algorithm.
 *
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class Run {
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_A = "hsa.faa";
	private static final String SAMPLE_PC_A = "hsa_pc.faa";
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_B = "cel.faa";
	private static final String SAMPLE_PC_B = "cel_pc.faa";
	
	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(Run.class.getName());
	private static final Options options = new Options();

	private static Options getOptions()
	{
		// general options
		options.addOption("h", false, "show this help message");
		options.addOption("v", true, "verbosity: 0 (default): print only warnings;\n1: print info messages");
		options.addOption("a", false, "print alignment");
		options.addOption("b", true, "bitscore cutoff");
		
		// debugging / negative control options
		options.addOption("D", false, "run debugging examples");
		options.addOption("A", false, "use adjusted BLOSUM62 matrix at coiled-coil positions");
		options.addOption("B", false, "use BLOSUM62 matrix at coiled-coil positions");
		options.addOption("LR", false, "load rounded matrix");
		options.addOption("R", false, "round matrix after scaling");
		options.addOption("P", false, "plain Smith-Waterman: do not use coiled-coil correction");
		options.addOption("N", false, "negative control: invert order of coiled-coil matrices");
		
		// db options
		options.addOption("p1", true, "protein sequences 1");
		options.addOption("c1", true, "coiled-coil prediction for protein sequences 1");
		options.addOption("s1", true, "optional regex for a sequence to use out of db1");
		options.addOption("p2", true, "protein sequences 1");
		options.addOption("c2", true, "coiled-coil prediction for protein sequences 2");
		options.addOption("s2", true, "optional regex for a sequence to use out of db2");
		options.addOption("s", false, "symmetric input: don't need p2/c2");

		// S-W params
		options.addOption("PE", true, "parameter: gap extension penalty");
		options.addOption("PO", true, "parameter: gap open penalty");
		options.addOption("PC", true, "parameter: coiled coil mismatch penalty");

		// parse file name to determine parameters
		options.addOption("F", true, "infer parameters from filename: keywords: blosum, pc0.1");
		
		return options;
	}
	
	private static void printUsage()
	{
		HelpFormatter h = new HelpFormatter();
		h.printHelp("ccaligner -p1 proteins1.fasta -c1 coil_predicton1.fasta -p2 proteins1.fasta -c2 coil_predicton1.fasta", options);
	}
	
	private static boolean checkDB(String pc_name, HashRichSequenceDB db, HashRichSequenceDB pcdb)
	{
    	RichSequenceIterator it = db.getRichSequenceIterator();
    	boolean error_found = false;
    	while (it.hasNext())
    	{
    		try
    		{
    			RichSequence s1 = it.nextRichSequence();
	    		String name = s1.getName(); 
	    		try
	    		{
	    			RichSequence s2 = pcdb.getRichSequence(name);
	    			if (s1.length() != s2.length())
	    			{
		    			System.err.println("different lengths in "+pc_name+": "+name);
		    			error_found = true; 
	    			}
	    		} catch (IllegalIDException e) {
	    			System.err.println("missing from "+pc_name+": "+name);
	    			error_found = true; 
	    		}
    		} catch (Exception e)
    		{
    			e.printStackTrace();
    			error_found = true;
    		}
    	}
    	return error_found;	
	}
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
        try {
        	logger.getParent().setLevel(Level.WARNING);
			DecimalFormat f1 = new DecimalFormat("0.00");

        	CommandLineParser parser = new GnuParser();
        	CommandLine cmd = parser.parse( getOptions(), args);
        	
        	if(cmd.hasOption("h"))
        	{
        		printUsage();
        		return;
        	}
        	
        	if(cmd.hasOption("v"))
        	{
        		String v = cmd.getOptionValue("v");
        		if (0 != v.compareTo("0"))
        		{
        			logger.getParent().setLevel(Level.INFO);
				}
        	}
        	
        	HashRichSequenceDB db1, db2, pcdb1, pcdb2;

        	boolean symm = cmd.hasOption("s");
        	
        	if(cmd.hasOption("D"))
        	{
            	logger.info("Running example...");
            	db1 = loadSequences(SAMPLE_SEQUENCE_A, cmd.getOptionValue("s1", ""));  
            	db2 = loadSequences(SAMPLE_SEQUENCE_B, cmd.getOptionValue("s2", ""));
            	pcdb1 = loadSequences(SAMPLE_PC_A, cmd.getOptionValue("s1", ""));
            	pcdb2 = loadSequences(SAMPLE_PC_B, cmd.getOptionValue("s2", ""));
        	}
        	else
        	{
        		// check if all needed params are there
        		ArrayList<String> missing = new ArrayList<String>();
        		for (String opt : "p1 p2 c1 c2".split(" "))
        		{
        			if (symm && opt.endsWith("2")) continue;
        			if (!cmd.hasOption(opt)) missing.add(opt);
        		}
        		
        		if (missing.size() > 0)
        		{
        			String ls = "";
        			for (String s : missing) ls += s + " ";
        			System.err.println("Please specify the following options for files to use: " + ls);
        			printUsage();
        			return;
        		}
        		
            	db1 = loadSequences(cmd.getOptionValue("p1"), cmd.getOptionValue("s1", ""));
            	pcdb1 = loadSequences(cmd.getOptionValue("c1"), cmd.getOptionValue("s1", ""));
            	
            	if (symm)
            	{
            		db2 = db1;
            		pcdb2 = pcdb1;
            	}
            	else
            	{
                	db2 = loadSequences(cmd.getOptionValue("p2"), cmd.getOptionValue("s2", ""));
                	pcdb2 = loadSequences(cmd.getOptionValue("c2"), cmd.getOptionValue("s2", ""));
            	}

            	// check dbs for consistency
            	boolean err1 = checkDB(cmd.getOptionValue("c1"),db1,pcdb1);
            	boolean err2 = symm ? err1 : checkDB(cmd.getOptionValue("c2"),db2,pcdb2);
            	if (err1 || err2) return;
        	} 
        	
        	float paramGapOpen = 10.0f;
        	if (cmd.hasOption("PO")) paramGapOpen = Float.valueOf(cmd.getOptionValue("PO")).floatValue();

        	float paramGapExt = 1.0f;
        	if (cmd.hasOption("PE")) paramGapExt = Float.valueOf(cmd.getOptionValue("PE")).floatValue();

        	float paramCoil = 0.15f;
        	if (cmd.hasOption("PC")) paramCoil = Float.valueOf(cmd.getOptionValue("PC")).floatValue();

        	ArrayList<Matrix> matrices = null;
        	
        	boolean coiled_coil_sw = !cmd.hasOption("P");
        	boolean zero_matrix = cmd.hasOption("N");
        	boolean blosum_matrix = cmd.hasOption("B");
        	boolean exact_matrix = !cmd.hasOption("LR");
        	boolean round_matrix = cmd.hasOption("R");
        	int adjusted_matrix = 0;
        	
        	float bitscore_cutoff = Float.valueOf(cmd.getOptionValue("b", "10"));
        	
        	if (cmd.hasOption("A")) adjusted_matrix = 1;
        	
        	if (cmd.hasOption("F"))
        	{
        		for (String token: cmd.getOptionValue("F").split("_"))
        		{
        			if (token.equals("ccalign")) { coiled_coil_sw = true; }
        			else if (token.equals("swalign")) { coiled_coil_sw = false; }
        			else if (token.equals("blosum")) { blosum_matrix = true; }
        			else if (token.equals("blosumx")) { blosum_matrix = true; exact_matrix = true; }
        			else if (token.equals("blosumxr")) { blosum_matrix = true; exact_matrix = true; round_matrix = true; }
        			else if (token.equals("ccx")) { exact_matrix = true; }
        			else if (token.equals("adjusted")) { adjusted_matrix = 1; }
        			else if (token.equals("adjustedx")) { adjusted_matrix = 1; exact_matrix = true; }
        			else if (token.equals("adjustedxr")) { adjusted_matrix = 1; exact_matrix = true; round_matrix = true; }
        			else if (token.equals("adjusted2")) { adjusted_matrix = 2; }
        			else if (token.equals("adjusted2x")) { adjusted_matrix = 2; exact_matrix = true; }
        			else if (token.equals("adjusted2xr")) { adjusted_matrix = 2; exact_matrix = true; round_matrix = true; }
        			else if (token.equals("zero")) { zero_matrix = true; }
        			else if (token.startsWith("pc")) { paramCoil = Float.valueOf(token.substring(2)).floatValue(); }
        			else 
        			{
        				throw new Exception("unknown token in command line: "+token);
        			}
        		}
        	}
        	
        	System.out.println("# coiled-coil mismatch penalty: " + paramCoil);
        	if (!coiled_coil_sw) { System.out.println("# plain Smith-Waterman search"); }
        	if (zero_matrix) { System.out.println("# using zero coiled-coil matrix"); }
        	if (blosum_matrix) { System.out.println("# using BLOSUM matrix"); }
        	if (adjusted_matrix > 0) { System.out.println("# using adjusted BLOSUM matrix"); }
        	if (exact_matrix) { System.out.println("# loading exact (not rounded) matrix"); }
        	if (round_matrix) { System.out.println("# rounding matrix"); }
        	System.out.println("# bitscore cutoff: " + f1.format(bitscore_cutoff));
        	
        	String blosum_fn = "BLOSUM62";
        	if (exact_matrix) blosum_fn += "x"; 
        	
        	Matrix blosum = MatrixLoader.load(blosum_fn);
        	if (exact_matrix) blosum.scaleScores(2);
    		if (round_matrix) blosum.roundScores();

    		// load coiled coil matrices unless we want to use plain S-W for control purposes
        	if (coiled_coil_sw)
        	{
        		matrices = new ArrayList<Matrix>();
        	
	        	for (char c : "abcdefg".toCharArray())
	        	{
	        		Matrix matrix;
	        		if (zero_matrix)
	        		{
	        			// use empty matrix
	        			matrix = new Matrix();
	        		}
	        		else if (blosum_matrix || adjusted_matrix > 0)
	        		{
	        			// load BLOSUM matrix, optionally use exact (instead of rounded) version
	        			matrix = MatrixLoader.load(blosum_fn);
	        			if (adjusted_matrix > 0)
	        			{
	        				float blosum_score = 0.6979f;
	        				float ad_score = (0.2466f + 0.2373f) / 2;
	        				float eg_score = (0.1848f + 0.1891f) / 2;
	        				float bcf_score = (0.0736f + 0.0926f + 0.0530f) / 3;
	        				
	        				float scale = 1.0f / ((adjusted_matrix == 1) ? ad_score : blosum_score); 
	        				if (c == 'a' || c == 'd') scale *= ad_score; 
	        				else if (c == 'e' || c == 'g') scale *= eg_score;
	        				else scale *= bcf_score;
	        				
	        				matrix.scaleScores(scale);
	        			}
	                	if (exact_matrix) matrix.scaleScores(2);
	        		} 
	        		else if (exact_matrix)
	        		{
	        			matrix = MatrixLoader.load(c + "_blosum.sij");
	        			matrix.scaleScores(2);
	        		}
	        		else
	        		{
	        			matrix = MatrixLoader.load(c + "_blosum.iij");
	        		}
	        		
	        		if (round_matrix) matrix.roundScores();
	        		
	        		matrices.add(matrix);
	        	}
        	}
        	
        	// set up estimates for remaining time
        	int sum1 = 0, sum2 = 0;

        	RichSequenceIterator it1 = db1.getRichSequenceIterator();
        	while (it1.hasNext())
        	{
        		RichSequence s1 = it1.nextRichSequence();
        		sum1 += s1.length();
        	}

        	BigInteger total_todo;
        	
        	it1 = db2.getRichSequenceIterator();
        	while (it1.hasNext())
        	{
        		RichSequence s1 = it1.nextRichSequence();
        		sum2 += s1.length();
        	}
        	total_todo = BigInteger.valueOf(sum1).multiply(BigInteger.valueOf(sum2)); 
	        	
        	BigInteger total_done = BigInteger.valueOf(0);
    		long start = System.currentTimeMillis();
    		long last_notification = start - 9000; // print first notification after 1 second 
    		System.err.println(total_todo);
    		// PE 1, PO 10
    		float lambda = 0.2288567f;
    		float kappa = 0.3101115f;
    		
    		// perform S-W alignments
        	it1 = db1.getRichSequenceIterator();
        	while (it1.hasNext())
        	{
        		RichSequence s1 = it1.nextRichSequence();
        		RichSequence pc1 = pcdb1.getRichSequence(s1.getName());

        		assert s1.length() == pc1.length() : s1.getName();
        		
    			RichSequenceIterator it2 = db2.getRichSequenceIterator();
        		while (it2.hasNext())
            	{
        			RichSequence s2 = it2.nextRichSequence();
        			
        			total_done = total_done.add(BigInteger.valueOf(s1.length()*s2.length()));

        			// in the symmetrical case, only do upper triangle
        			if (symm && (s1.getName().compareTo(s2.getName()) < 0)) continue;
        			
        			RichSequence pc2 = pcdb2.getRichSequence(s2.getName());
	
        			try
        			{
        				Alignment alignment = SmithWatermanGotoh.align(s1, s2, pc1, pc2, matrices, blosum, paramGapOpen, paramGapExt, paramCoil);

            			if (cmd.hasOption("a"))
            			{
    		    	        System.out.println ( alignment.getSummary() );
    		    	        System.out.println ( new Pair().format(alignment) );
    	
    		    	        System.out.println ( ">"+alignment.getName1() );
    		    	        System.out.println ( alignment.getSequence1() );
    		    	        System.out.println ( ">"+alignment.getName2() );
    		    	        System.out.println ( alignment.getSequence2() );
            			}
            			
            			float score = alignment.getScore();
            			double bitscore = (lambda * score - Math.log(kappa)) / Math.log(2);
            			// we cannot calculate a comparable e-value to blast as blast uses an "effective search space"
            			
            			if (bitscore >= bitscore_cutoff)
            			{
                			String result = s1.getName()+ "\t"+s2.getName()+"\t"+f1.format(bitscore)+"\t"+f1.format(100.0*alignment.getIdentity()/alignment.getSequence1().length)+"\t"
    							+ alignment.getStart1() + "\t" + (alignment.getStart1()+alignment.getSequence1().length)+"\t"
    							+ alignment.getStart2() + "\t" + (alignment.getStart2()+alignment.getSequence2().length);
                			System.out.println(result);
            			}
}
        			catch (Exception e)
        			{
        				System.err.println("Exception when aligning "+s1.getName()+" and "+s2.getName());
        				throw e;
        			}


	    	        
	    	        // print notification every 10 seconds on remaining time 
	        		long now = System.currentTimeMillis();
	        		
	        		if (now - last_notification > 10000)
	        		{
	        			long total_est = (total_todo.multiply(BigInteger.valueOf(now-start))).divide(total_done).longValue();
	        			long perc_done = BigInteger.valueOf(100).multiply(total_done).divide(total_todo).longValue();
	        			
	        			long millis = Math.round(total_est - (now-start));
	        			String left = String.format("%d:%02d:%02d",
	        					TimeUnit.MILLISECONDS.toHours(millis),
	        				    TimeUnit.MILLISECONDS.toMinutes(millis) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(millis)),
	        				    TimeUnit.MILLISECONDS.toSeconds(millis) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
	        				);
	        			
	        			System.err.println( perc_done + "% done, remaining: " + left);
	        			
	        			last_notification = now;
	        		}
            	}
        	}
        	
	        
	        logger.info("Finished running ccaligner");
	        
        } catch (Exception e) {
        	logger.log(Level.SEVERE, "Failed running ccaligner: " + e.getMessage(), e);
        	System.exit(1);
        }
    }
	
	/**
	 * 
	 * @param path location of the sequence
	 * @param filter 
	 * @return sequence string
	 * @throws IOException
	 * @throws BioException 
	 */
	private static HashRichSequenceDB loadSequences(String path, String filter) throws IOException, BioException {
		
		Reader reader = null;
		
		try
		{
			// try to read from file system
			reader = new FileReader(path);
		} catch (FileNotFoundException e)
		{
			// if the file doesn't exist, check if it is an internal example
			try
			{
				reader = new InputStreamReader(Run.class.getClassLoader().getResourceAsStream("ccaligner/run/sequences/"+path));
			} catch (Exception ee)
			{
				throw e; 
			}
		}
		
	    Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
	    SimpleNamespace ns = new SimpleNamespace("biojava");

	    HashRichSequenceDB db = new HashRichSequenceDB();
	    
	    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(new BufferedReader(reader), alpha.getTokenization("token"), ns);
	    while (iterator.hasNext()) {
	    	RichSequence seq = iterator.nextRichSequence();
	    	if (filter.isEmpty() || seq.getName().matches(filter)) db.addRichSequence(seq);
	    }
	    
		return db;
	}

}