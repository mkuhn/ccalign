package ccaligner.matrix;

public class MatrixCoord {

	private final int i;
	private final int j;

	public MatrixCoord(String coord) {
		this.i = coord.charAt(0);
		this.j = coord.charAt(1);
	}

	public int getI() {
		return i;
	}

	public int getJ() {
		return j;
	}
}
