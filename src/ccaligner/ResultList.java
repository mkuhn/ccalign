package ccaligner;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

public class ResultList implements Iterable<AlignmentResult> {

	private TreeSet<AlignmentResult> results;

	public ResultList()
	{
		results = new TreeSet<AlignmentResult>();
	}
	
	public Iterator<AlignmentResult> iterator() {
		return results.iterator();
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
	
}
