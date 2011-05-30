package ccaligner;

import java.text.DecimalFormat;

import ccaligner.util.Commons;

public class AlignmentResult implements Comparable<AlignmentResult> {

	private static final DecimalFormat f1 = new DecimalFormat("0.00");

	private final String p1;
	private final String p2;
	private final float bitscore;
	private final float identity;
	private final int start1;
	private final int end1;
	private final int start2;
	private final int end2;
	private String method;
	private final String message;
	
	public AlignmentResult(String p1, String p2, String message) {
		this.p1 = p1;
		this.p2 = p2;
		this.bitscore = this.identity = this.start1 = this.end1 = this.start2 = this.end2 = 0;
		this.method = null;
		this.message = message;
	}
	
	public AlignmentResult(String p1, String p2, float bitscore, float identity, int start1,
			int end1, int start2, int end2, String method) throws Exception {
		this.p1 = Commons.extractName(p1);
		this.p2 = Commons.extractName(p2);
		this.bitscore = bitscore;
		this.identity = identity;
		this.start1 = start1;
		this.end1 = end1;
		this.start2 = start2;
		this.end2 = end2;
		this.method = method;
		this.message = null;
	}
	
	public AlignmentResult(Alignment alignment)
	{
		this.p1 = alignment.getName1();
		this.p2 = alignment.getName2();
		this.bitscore = alignment.getBitscore();
		this.identity = (float) (100.0*alignment.getIdentity()/alignment.getSequence1().length);
		this.start1 = alignment.getStart1();
		this.end1 = alignment.getStart1()+alignment.getSequence1().length;
		this.start2 = alignment.getStart2();
		this.end2 = alignment.getStart2()+alignment.getSequence2().length;
		this.method = "CC";
		this.message = null;
	}
	
	public AlignmentResult(String line) throws Exception
	{
		String[] fields = line.split("\t");
		
		this.p1 = Commons.extractName(fields[0]);
		this.p2 = Commons.extractName(fields[1]);
		
		if (fields.length == 3)
		{
			this.bitscore = this.identity = this.start1 = this.end1 = this.start2 = this.end2 = 0;
			this.method = null;
			this.message = fields[2];			
		}
		else
		{
			this.bitscore = Float.valueOf(fields[2]);
			this.identity = Float.valueOf(fields[3]);
			this.start1 = Integer.valueOf(fields[4]);
			this.end1 = Integer.valueOf(fields[5]);
			this.start2 = Integer.valueOf(fields[6]);
			this.end2 = Integer.valueOf(fields[7]);
			this.method = (fields.length > 8) ? fields[8] : "SW";
			this.message = null;
		}
	}

	public int compareTo(AlignmentResult other) {
		// higher bitscore should be listed first
		int result = Float.compare(other.bitscore, bitscore);
		if (result != 0) return result;
		
		result = p1.compareTo(other.p1); 
		if (result != 0) return result;

		return p2.compareTo(other.p2);
	}

	public String toString()
	{
		if (message != null) return p1 + "\t"+ p2 +"\t" + message;
		
		return p1 + "\t"+ p2 +"\t"+f1.format(bitscore)+"\t"+f1.format(identity)+"\t"
				+ (1+start1) + "\t" + end1 + "\t" + (1+start2) + "\t" + end2 + "\t" + method;
	}
	
	public boolean isRecomputed()
	{
		return method.contentEquals("CC");
	}
	
	public void setMethod(String m)
	{
		method = m;
	}
	
	
	public float getBitscore()
	{
		return bitscore;
	}
	
	public String getMessage()
	{
		return message;
	}
	
	public String getName1()
	{
		return p1;
	}

	public String getName2()
	{
		return p2;
	}
}
