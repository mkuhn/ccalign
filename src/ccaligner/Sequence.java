package ccaligner;

public class Sequence {
	public final String name;
	public final Residue[] residues;

	public Sequence(String name, Residue[] residues) {
		this.name = name;
		this.residues = residues;
	}
}
