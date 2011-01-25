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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

import ccaligner.Alignment;
import ccaligner.AlignmentResult;
import ccaligner.Residue;
import ccaligner.ResultList;
import ccaligner.Sequence;
import ccaligner.SmithWatermanGotoh;
import ccaligner.formats.Pair;
import ccaligner.matrix.Matrix;
import ccaligner.matrix.MatrixLoader;
import ccaligner.util.Commons;

/**
 * Example of using JAligner API to align P53 human against
 * P53 mouse using Smith-Waterman-Gotoh algorithm.
 *
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

class DoRun // implements Runnable
{
	private Sequence seq1;
	private Sequence seq2;
	private float paramGapOpen;
	private float paramGapExt;
	private float paramCoilMatch;
	private float paramCoilMismatch;
	private ArrayList<Matrix> matrices;
	private Matrix blosum;
	private boolean print_alignment;
	
	public DoRun(Sequence seq1, Sequence seq2, float paramGapOpen,
			float paramGapExt, float paramCoilMatch, float paramCoilMismatch, ArrayList<Matrix> matrices,
			Matrix blosum, boolean print_alignment) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.paramGapOpen = paramGapOpen;
		this.paramGapExt = paramGapExt;
		this.paramCoilMatch = paramCoilMatch;
		this.paramCoilMismatch = paramCoilMismatch;
		this.matrices = matrices;
		this.blosum = blosum;
		this.print_alignment = print_alignment;
	}



	public AlignmentResult run() throws Exception
	{
		try
		{
			Alignment alignment = SmithWatermanGotoh.align(seq1, seq2, matrices, blosum, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch);

			if (print_alignment)
			{
    	        System.out.println ( alignment.getSummary() );
    	        System.out.println ( new Pair().format(alignment) );

    	        System.out.println ( ">"+alignment.getName1() );
    	        System.out.println ( alignment.getSequence1() );
    	        System.out.println ( ">"+alignment.getName2() );
    	        System.out.println ( alignment.getSequence2() );
			}

			return new AlignmentResult(alignment);
		}
		catch (OutOfMemoryError e)
		{
			// if the sequences are too big to compare with the current memory limits, take not of this 
			// for later re-calculation
			System.err.println("Not enough memory to align "+seq1.name+" and "+seq2.name);
			return new AlignmentResult(seq1.name, seq2.name, "ERROR: Out of memory");
		}
		catch (Exception e)
		{
			System.err.println("Exception when aligning "+seq1.name+" and "+seq2.name);
			throw e;
		}
	}
}

public class Run {
	
	private static final BigInteger big0 = BigInteger.valueOf(0);
	private static final BigInteger big1 = BigInteger.valueOf(1);
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_A = "hsa.faa";
	private static final String SAMPLE_PC_A = "hsa.tsv";
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_B = "cel.faa";
	private static final String SAMPLE_PC_B = "cel.tsv";
	
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
		options.addOption("r", true, "read previous (Smith-Waterman) results from this file (or stdin if the parameter is '--')");
		options.addOption("rn", true, "the number of top hits that should be recomputed (in conjunction with -r)");
		options.addOption("rp", true, "1 or 2: recompute first or second protein row, not complete matrixl;\n-1: compute scores for missing proteins, e.g. due to out-of-memory errors");
		
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
		options.addOption("PC", true, "parameter: coiled coil match reward");
		options.addOption("PX", true, "parameter: coiled coil mismatch penalty");

		// parse file name to determine parameters
		options.addOption("F", true, "infer parameters from filename: keywords: blosum, pc0.1");
		
		return options;
	}
	
	private static void printUsage()
	{
		HelpFormatter h = new HelpFormatter();
		h.printHelp("ccaligner -p1 proteins1.fasta -c1 coil_predicton1.fasta -p2 proteins1.fasta -c2 coil_predicton1.fasta", options);
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
        	
        	Map<String,Sequence> seqs1;
        	Map<String,Sequence> seqs2;

        	boolean symm = cmd.hasOption("s");
        	
        	if(cmd.hasOption("D"))
        	{
            	logger.info("Running example...");
            	seqs1 = loadSequences(SAMPLE_SEQUENCE_A, SAMPLE_PC_A, cmd.getOptionValue("s1", ""));  
            	seqs2 = loadSequences(SAMPLE_SEQUENCE_B, SAMPLE_PC_B, cmd.getOptionValue("s2", ""));
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
        		
            	seqs1 = loadSequences(cmd.getOptionValue("p1"), cmd.getOptionValue("c1"), cmd.getOptionValue("s1", ""));
            	
            	if (symm)
            	{
            		seqs2 = seqs1;
            	}
            	else
            	{
                	seqs2 = loadSequences(cmd.getOptionValue("p2"), cmd.getOptionValue("c2"), cmd.getOptionValue("s2", ""));
            	}
        	} 
        	
        	float paramGapOpen = 10.0f;
        	if (cmd.hasOption("PO")) paramGapOpen = Float.valueOf(cmd.getOptionValue("PO")).floatValue();

        	float paramGapExt = 1.0f;
        	if (cmd.hasOption("PE")) paramGapExt = Float.valueOf(cmd.getOptionValue("PE")).floatValue();

        	float paramCoilMatch = 0.15f;
        	if (cmd.hasOption("PC")) paramCoilMatch = Float.valueOf(cmd.getOptionValue("PC")).floatValue();

        	float paramCoilMismatch = paramCoilMatch;
        	if (cmd.hasOption("PX")) paramCoilMismatch = Float.valueOf(cmd.getOptionValue("PX")).floatValue();

        	ArrayList<Matrix> matrices = null;
        	
        	boolean coiled_coil_sw = !cmd.hasOption("P");
        	boolean zero_matrix = cmd.hasOption("N");
        	boolean blosum_matrix = cmd.hasOption("B");
        	boolean exact_matrix = !cmd.hasOption("LR");
        	boolean round_matrix = cmd.hasOption("R");
        	int adjusted_matrix = 0;
        	boolean print_alignment = cmd.hasOption("a");
        	
        	float bitscore_cutoff = Float.valueOf(cmd.getOptionValue("b", "10"));
        	
        	if (cmd.hasOption("A")) adjusted_matrix = 1;
        	
        	if (cmd.hasOption("F"))
        	{
        		boolean mismatch_set = false;
        		boolean match_set = false;	
        		
        		// parse individual options
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
        			else if (token.startsWith("pc")) { paramCoilMatch = Float.valueOf(token.substring(2)).floatValue(); match_set = true; }
        			else if (token.startsWith("px")) { paramCoilMismatch = Float.valueOf(token.substring(2)).floatValue(); mismatch_set = true;  }
        			else 
        			{
        				throw new Exception("unknown token in command line: "+token);
        			}
        		}
        		
        		// special case: if no mismatch penalty has been specified, use match reward
        		if (match_set && !mismatch_set)
        		{
        			paramCoilMismatch = paramCoilMatch;
        		}
        	}
        	
        	System.out.println("# coiled-coil match reward: " + paramCoilMatch);
        	System.out.println("# coiled-coil mismatch penalty: " + paramCoilMismatch);
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

        	for (Sequence s : seqs1.values())
        	{
        		sum1 += s.residues.length;
        	}

        	for (Sequence s : seqs2.values())
        	{
        		sum2 += s.residues.length;
        	}

        	// if this option is set, an existing set of scores is re-computed to avoid running
        	// CCAlign on too many non-relevant proteins
        	if (cmd.hasOption("r"))
    		{
            	BufferedReader br = new BufferedReader(openFile(cmd.getOptionValue("r")));
            	String line;

            	int to_check = 10;
            	if (cmd.hasOption("rn")) to_check = Integer.valueOf(cmd.getOptionValue("rn"));

            	int recompute_pass = 0;
            	if (cmd.hasOption("rp")) recompute_pass = Integer.valueOf(cmd.getOptionValue("rp"));

            	Map<String,ResultList> results1;
            	Map<String,ResultList> results2; 
            	
            	if (recompute_pass == 0)
            	{
            		// keep matrix in memory
            		results1 = new HashMap<String,ResultList>(seqs1.size());
            		results2 = new HashMap<String,ResultList>(seqs2.size());
            	}
            	else
            	{
            		// don't keep matrix in memory, still need this for temporary storage
            		results1 = new HashMap<String,ResultList>();
            		results2 = new HashMap<String,ResultList>();
            	}
            	
        		long start = System.currentTimeMillis();
        		long last_notification = start - 9000; // print first notification after 1 second 

            	BigInteger total_todo = BigInteger.valueOf(seqs1.size());
            	
            	if (recompute_pass == 2) total_todo = BigInteger.valueOf(seqs2.size());
            	
            	BigInteger total_done = big0;
        		ResultList rl = null;
        		
        		String last_name = "";
        		
            	while ((line = br.readLine()) != null)
            	{
            		// recompute_pass -1: just echo and recompute lines with errors
            		if (line.startsWith("#")) 
            		{
            			if (recompute_pass == -1) System.out.println(line); 
            			continue;
            		}
            		
            		AlignmentResult ar = new AlignmentResult(line);
            		
            		if (ar.getBitscore() < bitscore_cutoff) continue;
            		
            		// for recompute_pass 1 or 2, recompute for first or second column 
            		String name = (recompute_pass < 2) ? ar.getName1() :  ar.getName2();

        			if (recompute_pass > 0)
        			{
        				if (!name.contentEquals(last_name))
        				{
        					if (rl != null)
        					{
	        					results1.put("", rl);
	            				recompute(results1, null, recompute_pass, bitscore_cutoff, to_check, seqs1, seqs2, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, matrices, blosum, 0);
	            				results1.clear();
		            			total_done = total_done.add(big1);
        					}
                			rl = new ResultList();
                			last_name = name;
        				}

        				rl.add(ar);
        			}
        			else if (recompute_pass == -1)
        			{
        				// only re-compute missing lines
        				if (ar.getMessage() == null)
        				{
        					System.out.println(line);
        				}
        				else
        				{
                			DoRun task = new DoRun(seqs1.get(ar.getName1()), seqs2.get(ar.getName2()), paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, matrices, blosum, print_alignment);
                			AlignmentResult result = task.run();
                	        if (result.getBitscore() >= bitscore_cutoff) System.out.println(result.toString());
        				}
        			}
        			else
        			{
        				// re-compute everything at once
                		rl = results1.get(name);

	            		if (rl == null)
	            		{
	               			rl = new ResultList();
	               			results1.put(name, rl);
	            			total_done = total_done.add(big1);
	            		}

	            		rl.add(ar);
        			}
        			
            		
            		last_notification = printProgress(total_todo, total_done, last_notification, start);
            	}
            	
    			if (recompute_pass == 0)
    			{
    				System.err.println("starting first pass through alignments, no output expected yet");
    				recompute(results1, results2, recompute_pass, bitscore_cutoff, to_check, seqs1, seqs2, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, matrices, blosum, sum1);
    				System.err.println("starting second pass through alignments, printing alignments");
    				recompute(results2, null, recompute_pass, bitscore_cutoff, to_check, seqs1, seqs2, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, matrices, blosum, sum2);
    			}
    		}
    		else
    		{
    			BigInteger total_done = BigInteger.valueOf(0);
        		long start = System.currentTimeMillis();
        		long last_notification = start - 9000; // print first notification after 1 second 

        		BigInteger total_todo = BigInteger.valueOf(sum1).multiply(BigInteger.valueOf(sum2)); 

            	// perform S-W alignments
            	for (Sequence seq1 : seqs1.values())
            	{
            		for (Sequence seq2 : seqs2.values())
                	{
            			total_done = total_done.add(BigInteger.valueOf(seq1.residues.length*seq2.residues.length));

            			// in the symmetrical case, only do upper triangle
            			if (symm && (seq1.name.compareTo(seq2.name) < 0)) continue;
            			
            			DoRun task = new DoRun(seq1, seq2, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, matrices, blosum, print_alignment);
            			AlignmentResult result = task.run();
            	        if (result.getBitscore() >= bitscore_cutoff || result.getMessage() != null) System.out.println(result.toString());
    	    	        
    	    	        last_notification = printProgress(total_todo, total_done, last_notification, start);
                	}
            	}    			
    		}
        	
        	System.out.println("#DONE");
	        logger.info("Finished running ccaligner");
	        
        } catch (Exception e) {
        	logger.log(Level.SEVERE, "Failed running ccaligner: " + e.getMessage(), e);
        	System.exit(1);
        }
    }
	
	private static long printProgress(BigInteger total_todo, BigInteger total_done, long last_notification, long start)
	{
		// print notification every 10 seconds on remaining time 
		long now = System.currentTimeMillis();
		
		if (now - last_notification > 10000 && total_done.compareTo(big0) != 0)
		{
			if (total_todo.compareTo(big0) == 1)
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
			}
			else
			{
				long millis = Math.round(now-start);
				String left = String.format("%d:%02d:%02d",
						TimeUnit.MILLISECONDS.toHours(millis),
					    TimeUnit.MILLISECONDS.toMinutes(millis) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(millis)),
					    TimeUnit.MILLISECONDS.toSeconds(millis) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
					);
				
				System.err.println( "processed " + total_done.toString() + " sequences in " + left);
			}
			
			last_notification = now;
		}
		
		return last_notification;
	}
	
	private static void recompute(Map<String,ResultList> results1, Map<String,ResultList> results2, int recompute_pass, float bitscore_cutoff, int to_check, Map<String,Sequence> seqs1, Map<String,Sequence> seqs2, 
			float paramGapOpen, float paramGapExt, float paramCoilMatch, float paramCoilMismatch, ArrayList<Matrix> matrices,
			Matrix blosum, int total_sequence_length) throws Exception
	{
		BigInteger total_done = BigInteger.valueOf(0);
		long start = 0, last_notification = 0; 

		if (total_sequence_length > 0)
		{
			start = System.currentTimeMillis();
			last_notification = start - 9000; // print first notification after 1 second
		}

		BigInteger total_todo = BigInteger.valueOf(total_sequence_length); 
		
		for (Entry<String, ResultList> entry : results1.entrySet())
		{
			ResultList rl = entry.getValue();
			String entry_name = entry.getKey();
			
			Collection<AlignmentResult> to_recompute;
			
			while (!(to_recompute = rl.removeFromTop(to_check)).isEmpty())
			{
				for (AlignmentResult ar : to_recompute)
				{
					Sequence seq1 = seqs1.get(ar.getName1());
					if (seq1 == null) throw new Exception("Cannot find sequence in seqs1: " + ar.getName1());
					Sequence seq2 = seqs2.get(ar.getName2());
					if (seq2 == null) throw new Exception("Cannot find sequence in seqs2: " + ar.getName2());

					// if the scores of the input are the same for non-CC proteins, can skip re-computing them
//					if (seq1.min_pvalue > 0.1 & seq2.min_pvalue > 0.1)
//					{
//						ar.setMethod("SW-no-cc");
//						rl.add(ar);
//					}
//					else
					{
	        			DoRun task = new DoRun(seq1, seq2, paramGapOpen, paramGapExt, paramCoilMatch, 
								paramCoilMismatch, matrices, blosum, false);
	        			
	        			ar = task.run();
	        			if (ar.getBitscore() >= bitscore_cutoff) rl.add(ar);
					}
				}
			}
			
			for (AlignmentResult ar : rl)
			{
				if (results2 != null)
				{
	        		String name = ar.getName2();
	        		ResultList rl2 = results2.get(name);
	        		if (rl2 == null)
	        		{
	        			rl2 = new ResultList();
	        			results2.put(name, rl2);
	        		}
	        		rl2.add(ar);
				}
				else
				{
					System.out.println(ar.toString());
				}
			}
			
			rl.clear();
			
			if (total_sequence_length > 0)
			{
				total_done = total_done.add(BigInteger.valueOf( (results2 == null ? seqs2 : seqs1).get(entry_name).residues.length ));
		        last_notification = printProgress(total_todo, total_done, last_notification, start);
			}
		}
	}
	
	
	/**
	 * 
	 * @param path location of the sequence
	 * @param filter 
	 * @return sequence string
	 * @throws Exception 
	 */
	private static Map<String, Sequence> loadSequences(String aa_path, String cc_path, String filter) throws Exception {
		
		Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
	    SimpleNamespace ns = new SimpleNamespace("biojava");

	    // first, read sequence from FASTA so that we know which sequences to align
	    Map<String,Integer> sequence_lengths = new HashMap<String,Integer>(); 
	    
	    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(new BufferedReader(openFile(aa_path)), alpha.getTokenization("token"), ns);
	    while (iterator.hasNext()) {
	    	RichSequence seq = iterator.nextRichSequence();
	    	
	    	if (filter.isEmpty() || seq.getName().matches(filter)) {
	    		sequence_lengths.put(seq.getName(), seq.length());
	    	}
	    }

	    // second, read coil predictions
	    Map<String,Sequence> sequences = new HashMap<String,Sequence>(); 

	    BufferedReader br = new BufferedReader(openFile(cc_path));
	    
	    ArrayList<Residue> residues = null;
	    String name = null;
	    
	    for (String line = br.readLine(); line != null; line = br.readLine()) {
	    	if (line.startsWith(">"))
	    	{
	    		if (residues != null)
	    		{
	    			sequences.put( name, new Sequence(name, residues.toArray(new Residue[0])) );
	    		}
	    		
	    		name = Commons.extractName(line);
	    		
	    		if (sequence_lengths.containsKey(name))
	    		{
	    			residues = new ArrayList<Residue>(sequence_lengths.get(name));
	    		}
	    		else
	    		{
	    			residues = null;
	    		}
	    	}
	    	else if (residues != null)
	    	{
	    		String[] l = line.split("\t");
	    		char residue = l[0].charAt(0); 
	    		int register = l[1].charAt(0) - 'a';
	    		float pv = Float.valueOf(l[2]);
	    		
	    		if (pv > 0.1)
	    		{
	    			register = -1;
	    		}
	    		
	    		residues.add( new Residue(residue, register, pv) );
	    	}
	    }
	    
		if (residues != null)
		{
			sequences.put( name, new Sequence(name, residues.toArray(new Residue[0])) );
		}
		
		for (String seq_name : sequence_lengths.keySet())
		{
			if (sequences.containsKey(seq_name)) continue;
			throw new Exception("Missing coiled-coil prediction for '"+seq_name+"' when reading from '"+cc_path+"'!");
		}

		return sequences;
	}
	
	

	private static Reader openFile(String path) throws IOException {
		
		if (path.contentEquals("--"))
		{
			return new InputStreamReader(System.in);
		}
		
		try
		{
			// try to read from file system
			return new FileReader(path);
		} catch (FileNotFoundException e)
		{
			// if the file doesn't exist, check if it is an internal example
			try
			{
				return new InputStreamReader(Run.class.getClassLoader().getResourceAsStream("ccaligner/run/sequences/"+path));
			} catch (Exception ee)
			{
				throw e; 
			}
		}
	}
}