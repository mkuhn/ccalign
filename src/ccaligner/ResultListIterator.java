/**
 * 
 */
package ccaligner;

import java.io.BufferedReader;
import java.util.Iterator;
import java.util.logging.Logger;

import ccaligner.matrix.Matrix;
import ccaligner.matrix.MatrixLoader;

/**
 * Given a result list file in tsv format, optionally with embedded matrices, read hits
 * and return them grouped by query.
 * 
 * @author mkuhn
 *
 */
public class ResultListIterator implements Iterator<ResultList> {

	
	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(ResultListIterator.class.getName());
	
	private BufferedReader br;
	private float bitscore_cutoff;
	private int recompute_pass;
	
	private boolean next_rl_valid = false;
	private ResultList next_rl = null;
	
	private AlignmentResult ar_cache = null;
	private String line_cache = null;
	
	private Exception exception = null;
	
	private MatrixLoader ml;
	
	private ResultList readNext() throws Exception
	{
		ResultList rl = new ResultList();
		
		String last_name = null;
		
    	while (true)
    	{
    		AlignmentResult ar = null;

    		if (ar_cache == null)
    		{
    			// no cached alignment: read new line from input, but check for cached line
    			String line;
    			if (line_cache == null) {
    				line = br.readLine();
    			} else {
    				line = line_cache;
    				line_cache = null;
    			}
        		
        		if (line == null)
        		{
        			return rl;
        		}

        		if (line.startsWith("#")) 
        		{
        			// Start of new comment section: return results so far
        			if (rl.size() > 0)
        			{
        				logger.fine("Returning new result list");
        				line_cache = line;
        				return rl;
        			}
        			
        			if (line.startsWith("## query"))
        			{
        				ml = new MatrixLoader(line);
        			}
        			else
        			{
        				if (ml == null)
        				{
            				logger.warning("Encountering adjusted matrix before matrix header, using standard header");
            				ml = new MatrixLoader();
        				}
        				Matrix m = ml.loadFromLine(line);
        				if (m != null)
        				{
        					rl.addMatrix(m.getId(), m);
        				}
            			// recompute_pass -1: just echo and recompute lines with errors
        				else if (recompute_pass == -1) System.out.println(line);
        			}
        			continue;
        		}

        		ar = new AlignmentResult(line);
        		if (ar.getBitscore() < bitscore_cutoff) continue;

    		}
    		else
    		{
    			// use left-over alignment result from last call
    			ar = ar_cache;
    			ar_cache = null;
    		}

    		String name = (recompute_pass < 2) ? ar.getName1() : ar.getName2();
    		if (last_name == null) last_name = name;
    		
			if (!name.contentEquals(last_name))
			{
				if (rl.size() > 0)
				{
					ar_cache = ar;
					return rl;
				}
				last_name = name;
			}

			if (ar != null) rl.add(ar);
    	}
	}
	
	
	public ResultListIterator(BufferedReader br, float bitscore_cutoff, int recompute_pass)
	{
		this.br = br;
		this.bitscore_cutoff = bitscore_cutoff;
		this.recompute_pass = recompute_pass;
	}
	
	public boolean hasNext() {
		if (!next_rl_valid)
			try {
				next_rl = readNext();
				next_rl_valid = true;
			} catch (Exception e) {
				exception = e;
				return false;
			}
		
		return next_rl.size() > 0;
	}

	public ResultList next() {
		if (!hasNext()) return null;
		next_rl_valid = false;
		return next_rl;
	}

	/**
	 * Optional method: not implemented
	 */
	public void remove() {}
	
	
	public Exception lastException()
	{
		return exception;
	}
}
