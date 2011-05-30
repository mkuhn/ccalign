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
import java.util.BitSet;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

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
import ccaligner.ResultListIterator;
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
	private final Sequence seq1;
	private final Sequence seq2;
	private final float paramGapOpen;
	private final float paramGapExt;
	private final float paramCoilMatch;
	private final float paramCoilMismatch;
	private final Matrix cc_matrix;
	private final Matrix mx_matrix;
	private final Matrix no_matrix;
	private final int adjusted_matrix;
	private final boolean print_alignment;
	
	public DoRun(Sequence seq1, Sequence seq2, float paramGapOpen,
			float paramGapExt, float paramCoilMatch, float paramCoilMismatch, Matrix cc_matrix, Matrix mx_matrix, Matrix no_matrix, int adjusted_matrix, boolean print_alignment) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.paramGapOpen = paramGapOpen;
		this.paramGapExt = paramGapExt;
		this.paramCoilMatch = paramCoilMatch;
		this.paramCoilMismatch = paramCoilMismatch;
		this.cc_matrix = cc_matrix;
		this.mx_matrix = mx_matrix;
		this.no_matrix = no_matrix;
		this.adjusted_matrix = adjusted_matrix;
		this.print_alignment = print_alignment;
	}



	public AlignmentResult run() throws Exception
	{
		try
		{
			Alignment alignment = SmithWatermanGotoh.align(seq1, seq2, cc_matrix, mx_matrix, no_matrix, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, adjusted_matrix);

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
	
	private static final DecimalFormat f1 = new DecimalFormat("0.00");

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

	private static final String SAMPLE_MATRICES_AB = "hsa_cel_matrix.tsv";

	private static Matrix blosum;
	
	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(Run.class.getName());
	private static final Options options = new Options();

	/**
	 * Options
	 */
	private static float paramGapOpen = 10.0f;
	private static float paramGapExt = 1.0f;
	private static float paramCoilMatch = 0.15f;
	private static float paramCoilMismatch = paramCoilMatch;
	private static boolean cc_comp_adj = true;
	private static boolean blosum_matrix = false;
	private static int adjusted_matrix = 1;
	private static boolean print_alignment = false;
	private static float bitscore_cutoff = 0;
	
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
		options.addOption("rx", false, "print warning for missing sequences (if not set: abort with error)(");
		
		// debugging / negative control options
		options.addOption("A", true, "use adjusted matrix at coiled-coil positions: 0 - no adjustment, 1 - group registers (ad, bcf, eg), 2 - individual registers");
		options.addOption("B", false, "use BLOSUM62 matrix at coiled-coil positions");
		options.addOption("C", true, "split protein into cc/non-cc regions for compositional matrix adjustment");
		options.addOption("D", false, "run debugging examples");
		
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
	
	private static Matrix getMatrix(String prefix, String name1, String name2, Map<String, Matrix> matrices, Matrix blosum)
	{
		if (!blosum_matrix)
		{
			if (cc_comp_adj)
			{
				if (prefix != "") prefix = prefix + "-";
				Matrix matrix = matrices.get(prefix+name1+"-"+prefix+name2);
				if (matrix != null) return matrix;
			}
	
			Matrix matrix = matrices.get(name1+"-"+name2); 
			if (matrix != null) return matrix; 
		}
		return blosum;
	}

	private static Matrix getMatrix(String prefix, String name1, String name2, ResultList rl, Matrix blosum)
	{
		return rl.getMatrix(prefix, name1, name2, blosum, blosum_matrix, cc_comp_adj);
	}

	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
        try {
        	logger.getParent().setLevel(Level.WARNING);
//        	logger.getParent().setLevel(Level.INFO);

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
        	
        	if (cmd.hasOption("PO")) paramGapOpen = Float.valueOf(cmd.getOptionValue("PO")).floatValue();

        	if (cmd.hasOption("PE")) paramGapExt = Float.valueOf(cmd.getOptionValue("PE")).floatValue();

        	if (cmd.hasOption("PC")) paramCoilMatch = Float.valueOf(cmd.getOptionValue("PC")).floatValue();

        	if (cmd.hasOption("PX")) paramCoilMismatch = Float.valueOf(cmd.getOptionValue("PX")).floatValue();

        	cc_comp_adj = Integer.valueOf(cmd.getOptionValue("C")) != 0;
        	blosum_matrix = cmd.hasOption("B");
        	print_alignment = cmd.hasOption("a");
        	
        	bitscore_cutoff = Float.valueOf(cmd.getOptionValue("b", "10"));
        	
        	if (cmd.hasOption("A")) adjusted_matrix = 1;
        	
        	if (cmd.hasOption("F"))
        	{
        		boolean mismatch_set = false;
        		boolean match_set = false;	
        		
        		// parse individual options
        		for (String token: cmd.getOptionValue("F").split("_"))
        		{
        			if (token.equals("nosplit")) { cc_comp_adj = false; }
        			else if (token.equals("blosum")) { blosum_matrix = true; }
        			else if (token.equals("adjusted")) { adjusted_matrix = 1; }
        			else if (token.equals("adjusted2")) { adjusted_matrix = 2; }
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
        	if (blosum_matrix) 
        	{ 
        		System.out.println("# using BLOSUM matrix"); 
        	}
        	else
        	{
        		if (cc_comp_adj) { System.out.println("# using compositional matrix adjustment for CC/non-CC parts"); }
        		else { System.out.println("# using compositional matrix adjustment for whole protein"); }
        	}
        	if (adjusted_matrix > 0) { System.out.println("# using adjusted BLOSUM matrix"); }
        	System.out.println("# bitscore cutoff: " + f1.format(bitscore_cutoff));
        	
        	blosum = MatrixLoader.load("BLOSUM62");
        	
    		Map<String,Matrix> matrices = new HashMap<String,Matrix>();
    		
    		if(cmd.hasOption("D"))
        	{
            	logger.info("Running example...");
            	seqs1 = loadSequences(SAMPLE_SEQUENCE_A, SAMPLE_PC_A, cmd.getOptionValue("s1", ""));  
            	seqs2 = loadSequences(SAMPLE_SEQUENCE_B, SAMPLE_PC_B, cmd.getOptionValue("s2", ""));
            	
            	matrices = MatrixLoader.loadMatrices(SAMPLE_MATRICES_AB);
            	
        	}
        	else
        	{
        		// check if all needed parameters are there
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
            	
            	int to_check = 10;
            	if (cmd.hasOption("rn")) to_check = Integer.valueOf(cmd.getOptionValue("rn"));

            	int recompute_pass = 0;
            	if (cmd.hasOption("rp")) recompute_pass = Integer.valueOf(cmd.getOptionValue("rp"));

            	boolean skip_missing = cmd.hasOption("rx");
            	
        		long start = System.currentTimeMillis();
        		long last_notification = start - 9000; // print first notification after 1 second 

            	BigInteger total_todo = BigInteger.valueOf(seqs1.size());
            	
            	if (recompute_pass == 2) total_todo = BigInteger.valueOf(seqs2.size());

            	ResultListIterator rli = new ResultListIterator(br, bitscore_cutoff, recompute_pass);
            	BigInteger total_done = big0;
            	
            	while (rli.hasNext())
            	{
            		ResultList rl = rli.next();
            		
            		if (recompute_pass == -1)
        			{
            			for (AlignmentResult ar : rl)
            			{
	        				// only re-compute missing lines
	        				if (ar.getMessage() == null)
	        				{
	        					System.out.println(ar.toString());
	        				}
	        				else
	        				{
	        					String name1 = ar.getName1();
	        					String name2 = ar.getName2();

	        					Matrix cc_matrix = getMatrix("cc", name1, name2, rl, blosum);
	        					Matrix mx_matrix = getMatrix("", name1, name2, rl, blosum);
	        					Matrix no_matrix = getMatrix("no", name1, name2, rl, blosum);
	
	                			DoRun task = new DoRun(seqs1.get(ar.getName1()), seqs2.get(ar.getName2()), paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, cc_matrix, mx_matrix, no_matrix, adjusted_matrix, print_alignment);
	                			AlignmentResult result = task.run();
	                	        if (result.getBitscore() >= bitscore_cutoff) System.out.println(result.toString());
	        				}
            			}
        			}
            		else
            		{
        				recompute(rl, recompute_pass, bitscore_cutoff, to_check, seqs1, seqs2, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, matrices, skip_missing);
            			total_done = total_done.add(big1);
            		}

            		last_notification = printProgress(total_todo, total_done, last_notification, start);
            	}
            	
            	if (rli.lastException() != null) throw rli.lastException();
    		}
    		else
    		{
    			// Run CCAlign for all pairs of query/subject
    			BigInteger total_done = BigInteger.valueOf(0);
        		long start = System.currentTimeMillis();
        		long last_notification = start - 9000; // print first notification after 1 second 

        		BigInteger total_todo = BigInteger.valueOf(sum1).multiply(BigInteger.valueOf(sum2)); 

            	// perform S-W alignments
            	for (Sequence seq1 : seqs1.values())
            	{
            		final String name1 = seq1.name;
            		for (Sequence seq2 : seqs2.values())
                	{
            			total_done = total_done.add(BigInteger.valueOf(seq1.residues.length*seq2.residues.length));

                		final String name2 = seq2.name;

            			// in the symmetrical case, only do upper triangle
            			if (symm && (name1.compareTo(name2) < 0)) continue;

    					Matrix cc_matrix = getMatrix("cc", name1, name2, matrices, blosum);
    					Matrix mx_matrix = getMatrix("", name1, name2, matrices, blosum);
    					Matrix no_matrix = getMatrix("no", name1, name2, matrices, blosum);

            			DoRun task = new DoRun(seq1, seq2, paramGapOpen, paramGapExt, paramCoilMatch, paramCoilMismatch, cc_matrix, mx_matrix, no_matrix, adjusted_matrix, print_alignment);
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
	
	private static void recompute(ResultList rl, int recompute_pass, float bitscore_cutoff, int to_check, Map<String,Sequence> seqs1, Map<String,Sequence> seqs2, 
			float paramGapOpen, float paramGapExt, float paramCoilMatch, float paramCoilMismatch, Map<String,Matrix> matrices,
			boolean skip_missing) throws Exception
	{
		Collection<AlignmentResult> to_recompute;
		
		while (!(to_recompute = rl.removeFromTop(to_check)).isEmpty())
		{
			for (AlignmentResult ar : to_recompute)
			{
				Sequence seq1 = seqs1.get(ar.getName1());
				if (seq1 == null)
				{
					if (skip_missing)
					{
						logger.warning("Missing sequence from seq1: " + ar.getName1());
						continue;
					}
					else
					{
						throw new Exception("Cannot find sequence in seqs1: " + ar.getName1());
					}
				}
				Sequence seq2 = seqs2.get(ar.getName2());
				if (seq2 == null)
				{
					if (skip_missing)
					{
						logger.warning("Missing sequence from seq2: " + ar.getName2());
						continue;
					}
					else
					{
						throw new Exception("Cannot find sequence in seqs2: " + ar.getName2());
					}
				}

				// if the scores of the input are the same for non-CC proteins, can skip re-computing them
				if (seq1.max_prob < 0.9 & seq2.max_prob < 0.9)
				{
					ar.setMethod("SW-no-cc");
					rl.add(ar);
				}
				else
				{
					String name1 = ar.getName1();
					String name2 = ar.getName2();
					Matrix cc_matrix = getMatrix("cc", name1, name2, rl, blosum);
					Matrix mx_matrix = getMatrix("", name1, name2, rl, blosum);
					Matrix no_matrix = getMatrix("no", name1, name2, rl, blosum);

					logger.info("Using CC Matrix: " + cc_matrix.getId());
					
        			DoRun task = new DoRun(seq1, seq2, paramGapOpen, paramGapExt, paramCoilMatch, 
							paramCoilMismatch, cc_matrix, mx_matrix, no_matrix, adjusted_matrix, false);
        			
        			ar = task.run();
        			if (ar.getBitscore() >= bitscore_cutoff) rl.add(ar);
				}
			}
		}
			
		for (AlignmentResult ar : rl)
		{
			System.out.println(ar.toString());
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
	    // note: we don't store the sequences here, but only the name and length
	    Map<String,Integer> sequence_lengths = new HashMap<String,Integer>(); 
	    
	    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(new BufferedReader(openFile(aa_path)), alpha.getTokenization("token"), ns);
	    while (iterator.hasNext()) {
	    	RichSequence seq = iterator.nextRichSequence();
	    	
	    	if (filter.isEmpty() || seq.getName().matches(filter)) {
	    		sequence_lengths.put(seq.getName(), seq.length());
	    	}
	    }

	    // second, read coil predictions and protein sequence from the coiled coil prediction
	    Map<String,Sequence> sequences = new HashMap<String,Sequence>(); 

	    BufferedReader br = new BufferedReader(openFile(cc_path));
	    
	    ArrayList<Residue> residues = null;
	    String name = null;
	    
	    for (String line = br.readLine(); line != null; line = br.readLine()) {
	    	// check here if we're interested in the prediction for this sequence
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
	    		
	    		if (l.length != 10)
				{
					throw new Exception("Unrecognized format of coiled-coil prediction for '"+name+"' in '"+cc_path+"':\n"+line);
				}
	    		
	    		char residue = l[0].charAt(0); 
	    		int register = l[1].charAt(0) - 'a';
	    		float prob = Float.valueOf(l[2]);
	    		
	    		if (prob < 0.9)
	    		{
	    			register = -1;
	    		}
	    		
	    		BitSet possible_registers = new BitSet(7);
	    		
	    		for (int i = 0; i < 7; i++)
	    			if (Float.valueOf(l[3+i]) >= 0.01) possible_registers.set(i); 
	    		
	    		residues.add( new Residue(residue, register, prob, possible_registers) );
	    	}
	    }
	    
		if (residues != null)
		{
			if (residues.size() != sequence_lengths.get(name))
			{
				throw new Exception("Coiled-coil prediction for '"+name+"' does not match size between '"+aa_path+"' and '"+cc_path+"'! ("+residues.size()+" vs. "+sequence_lengths.get(name)+")");
			}
			
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