package ccaligner.matrix;

import ccaligner.Residue;

public class MatrixCoord {

	private final short i;
	private final short j;

	public MatrixCoord(String coord) {
		this.i = Residue.aa_to_code[ coord.charAt(0) ];
		this.j = Residue.aa_to_code[ coord.charAt(1) ];
	}

	public short getI() {
		return i;
	}

	public short getJ() {
		return j;
	}
}
