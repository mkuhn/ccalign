package ccaligner;

public class Sequence {
	public final String name;
	public final Residue[] residues;
	public final float max_prob;

	public Sequence(String name, Residue[] residues) {
		this.name = name;
		this.residues = residues;
		
		float prob = 0;
		for (Residue r : residues)
		{
			if (r.cc_prob > prob) prob = r.cc_prob; 
		}
		max_prob = prob;
	}
}
