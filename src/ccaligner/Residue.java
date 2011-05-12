package ccaligner;

import java.util.BitSet;

public class Residue {

	public final char aa;
	public final int register;
	public final float cc_prob;
	public final BitSet possible_registers;

	public Residue(char aa, int register, float cc_prob, BitSet possible_registers) {
		this.aa = aa;
		this.register = register;
		this.cc_prob = cc_prob;
		this.possible_registers = possible_registers;
	}
	
}
