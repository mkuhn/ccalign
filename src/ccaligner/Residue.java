package ccaligner;

import java.util.BitSet;

public class Residue {

	public final char aa;
	public final short aa_code;
	public final int register;
	public final float cc_prob;
	public final BitSet possible_registers;

	public static final short[] aa_to_code;
	public static final char[] residues = "X-ABCDEFGHIKLMNPQRSTVWYZU*OJ".toCharArray();;
	
	static {
		 aa_to_code = new short[127];
		 
		 for (short i = 0; i < residues.length; i++)
			 aa_to_code[ residues[i] ] = i;
	}
	
	
	public Residue(char aa, int register, float cc_prob, BitSet possible_registers) {
		this.aa = aa;
		this.aa_code = aa_to_code[aa];
		this.register = register;
		this.cc_prob = cc_prob;
		this.possible_registers = possible_registers;
	}
	
}
