package ccaligner;

public class Residue {

	public final char aa;
	public final int register;
	public final float cc_pvalue;

	public Residue(char aa, int register, float cc_pvalue) {
		this.aa = aa;
		this.register = register;
		this.cc_pvalue = cc_pvalue;
	}
	
}
