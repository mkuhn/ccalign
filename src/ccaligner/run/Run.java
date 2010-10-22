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
import java.util.ArrayList;
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
		
		options.addOption("h", false, "show this help message");
		options.addOption("v", true, "verbosity: 0 (default): print only warnings;\n1: print info messages");
		options.addOption("D", false, "run debugging examples");
		options.addOption("p1", true, "protein sequences 1");
		options.addOption("c1", true, "coiled-coil prediction for protein sequences 1");
		options.addOption("p2", true, "protein sequences 1");
		options.addOption("c2", true, "coiled-coil prediction for protein sequences 2");

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
	    		String name = it.nextRichSequence().getName(); 
	    		try
	    		{
	    			pcdb.getRichSequence(name);
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

        	if(cmd.hasOption("D"))
        	{
            	logger.info("Running example...");
            	db1 = loadSequences(SAMPLE_SEQUENCE_A);  
            	db2 = loadSequences(SAMPLE_SEQUENCE_B);
            	pcdb1 = loadSequences(SAMPLE_PC_A);
            	pcdb2 = loadSequences(SAMPLE_PC_B);
        	}
        	else
        	{
        		ArrayList<String> missing = new ArrayList<String>();
        		for (String opt : "p1 p2 c1 c2".split(" "))
        		{
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
        		
            	db1 = loadSequences(cmd.getOptionValue("p1"));  
            	db2 = loadSequences(cmd.getOptionValue("p2"));
            	pcdb1 = loadSequences(cmd.getOptionValue("c1"));
            	pcdb2 = loadSequences(cmd.getOptionValue("c2"));

            	// check dbs for consistency
            	boolean err1 = checkDB(cmd.getOptionValue("c1"),db1,pcdb1);
            	boolean err2 = checkDB(cmd.getOptionValue("c2"),db2,pcdb2);
            	if (err1 || err2) return;
        	} 
        	
        	
        	ArrayList<Matrix> matrices = new ArrayList<Matrix>();
        	
        	for (char c : "abcdefg".toCharArray())
        	{
        		Matrix matrix = MatrixLoader.load(c + "_blosum.iij");
        		matrices.add(matrix);
        	}
        	
        	Matrix blosum = MatrixLoader.load("BLOSUM62");


        	ArrayList<String> results = new ArrayList<String>();
        	
        	RichSequenceIterator it1 = db1.getRichSequenceIterator();
        	while (it1.hasNext())
        	{
        		RichSequence s1 = it1.nextRichSequence();
        		RichSequence pc1 = pcdb1.getRichSequence(s1.getName());

        		assert s1.length() == pc1.length() : s1.getName();
        		
    			RichSequenceIterator it2 = db2.getRichSequenceIterator();
        		while (it2.hasNext())
            	{
        			RichSequence s2 = it2.nextRichSequence();
        			RichSequence pc2 = pcdb2.getRichSequence(s2.getName());
	
        			Alignment alignment = SmithWatermanGotoh.align(s1, s2, pc1, pc2, matrices, blosum, 10f, 0.5f);
	    	        
	    	        System.out.println ( alignment.getSummary() );
	    	        System.out.println ( new Pair().format(alignment) );

	    	        System.out.println ( ">"+alignment.getName1() );
	    	        System.out.println ( alignment.getSequence1() );
	    	        System.out.println ( ">"+alignment.getName2() );
	    	        System.out.println ( alignment.getSequence2() );

	    	        results.add(s1.getName()+"\t"+s2.getName()+"\t"+alignment.getScore()+"\t"+alignment.getCoilMatches());
            	}
        	}
        	
	        
	        logger.info("Finished running example");
	        
	        for (String s:results)
	        {
	        	System.out.println(s);
	        }
	        
        } catch (Exception e) {
        	logger.log(Level.SEVERE, "Failed running example: " + e.getMessage(), e);
        }
    }
	
	/**
	 * 
	 * @param path location of the sequence
	 * @return sequence string
	 * @throws IOException
	 * @throws BioException 
	 */
	private static HashRichSequenceDB loadSequences(String path) throws IOException, BioException {
		
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
	        db.addRichSequence(seq);
	    }
	    
		return db;
	}

}