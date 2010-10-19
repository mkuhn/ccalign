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
import jaligner.matrix.MatrixLoader;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.db.HashRichSequenceDB;
import org.biojavax.bio.seq.*;

/**
 * Example of using JAligner API to align P53 human aganist
 * P53 mouse using Smith-Waterman-Gotoh algorithm.
 *
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class Example {
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_A = "src/jaligner/example/sequences/asl.fasta";
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_B = "src/jaligner/example/sequences/t07c4.10.fasta";
	
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
	        
        	RichSequenceIterator it1 = db1.getRichSequenceIterator();
        	while (it1.hasNext())
        	{
        		RichSequence s1 = it1.nextRichSequence();

    			RichSequenceIterator it2 = db2.getRichSequenceIterator();
        		while (it2.hasNext())
            	{
        			RichSequence s2 = it2.nextRichSequence();
            		
        			Alignment alignment = SmithWatermanGotoh.align(s1, s2, MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
	    	        
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
	    
	    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
	            alpha.getTokenization("token"), ns);
	    while (iterator.hasNext()) {
	        db.addRichSequence(iterator.nextRichSequence());
	    }
	    
		return db;
	}

}