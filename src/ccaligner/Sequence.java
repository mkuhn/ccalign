package ccaligner;

public class Sequence {
	public final String name;
	public final Residue[] residues;
	public final float min_pvalue;

	public Sequence(String name, Residue[] residues) {
		this.name = name;
		this.residues = residues;
		
		float pv = 1;
		for (Residue r : residues)
		{
			if (r.cc_pvalue < pv) pv = r.cc_pvalue; 
		}
		min_pvalue = pv;
	}
}
