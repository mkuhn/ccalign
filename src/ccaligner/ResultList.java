package ccaligner;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import ccaligner.matrix.Matrix;

public class ResultList implements Iterable<AlignmentResult> {

	private TreeSet<AlignmentResult> results;
	private Map<String, Matrix> matrices;

	public ResultList()
	{
		results = new TreeSet<AlignmentResult>();
		matrices = new HashMap<String, Matrix>();
	}
	
	public Iterator<AlignmentResult> iterator() {
		return results.iterator();
	}

	public void addMatrix(String name, Matrix m)
	{
		matrices.put(name, m);
	}

	public boolean add(AlignmentResult r)
	{
		return results.add(r);
	}
	
	public Collection<AlignmentResult> removeFromTop(int to_check)
	{
		ArrayList<AlignmentResult> l = new ArrayList<AlignmentResult>(to_check);
		
		int checked = 0;
		
		for (AlignmentResult r : results)
		{
			if (!r.isRecomputed()) 
			{
				l.add(r);
			}
			if (++checked >= to_check) break; 
		}
		
		results.removeAll(l);
		
		return l;
	}
	
	public void clear()
	{
		results.clear();
	}
	
	public int size()
	{
		return results.size();
	}

	public Matrix getMatrix(String prefix, String name1, String name2, Matrix blosum, boolean blosum_matrix, boolean cc_comp_adj) 
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
}
