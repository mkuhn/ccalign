/*
 * $Id: Example.java,v 1.3 2005/04/03 19:38:21 ahmed Exp $
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

package jaligner.example;

import jaligner.Alignment;
import jaligner.SmithWatermanGotoh;
import jaligner.formats.Pair;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;

import java.io.*;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.db.HashRichSequenceDB;
import org.biojavax.bio.seq.*;

/**
 * Example of using JAligner API to align P53 human against
 * P53 mouse using Smith-Waterman-Gotoh algorithm.
 *
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class Example {
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_A = "src/jaligner/example/sequences/sass6.faa";
	private static final String SAMPLE_PC_A = "src/jaligner/example/sequences/sass6_pc.faa";
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_B = "src/jaligner/example/sequences/cel.faa";
	private static final String SAMPLE_PC_B = "src/jaligner/example/sequences/cel_pc.faa";
	
	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(Example.class.getName());
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
        try {
        	logger.info("Running example...");
        	
        	HashRichSequenceDB db1 = loadSequences(SAMPLE_SEQUENCE_A);  
        	HashRichSequenceDB db2 = loadSequences(SAMPLE_SEQUENCE_B);
        	HashRichSequenceDB pcdb1 = loadSequences(SAMPLE_PC_A);
        	HashRichSequenceDB pcdb2 = loadSequences(SAMPLE_PC_B);

        	ArrayList<Matrix> matrices = new ArrayList<Matrix>();
        	
        	for (char c : "abcdefg".toCharArray())
        	{
        		Matrix matrix = MatrixLoader.load(c + "_blosum.iij");
        		matrices.add(matrix);
        	}
        	
        	
        	Matrix blosum = MatrixLoader.load("BLOSUM62");
        	
        	RichSequenceIterator it1 = db1.getRichSequenceIterator();
        	while (it1.hasNext())
        	{
        		RichSequence s1 = it1.nextRichSequence();
        		RichSequence pc1 = pcdb1.getRichSequence(s1.getName());

    			RichSequenceIterator it2 = db2.getRichSequenceIterator();
        		while (it2.hasNext())
            	{
        			RichSequence s2 = it2.nextRichSequence();
        			RichSequence pc2 = pcdb2.getRichSequence(s2.getName());
	
        			Alignment alignment = SmithWatermanGotoh.align(s1, s2, pc1, pc2, matrices, blosum, 10f, 0.5f);
	    	        
	    	        System.out.println ( alignment.getSummary() );
	    	        System.out.println ( new Pair().format(alignment) );
            	}
        	}
        	
	        
	        logger.info("Finished running example");
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
	    BufferedReader br = new BufferedReader(new FileReader(path));
	    Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
	    SimpleNamespace ns = new SimpleNamespace("biojava");

	    HashRichSequenceDB db = new HashRichSequenceDB();
	    
	    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br, alpha.getTokenization("token"), ns);
	    while (iterator.hasNext()) {
	    	RichSequence seq = iterator.nextRichSequence();
	        db.addRichSequence(seq);
	    }
	    
		return db;
	}

}