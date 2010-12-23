/*
 * $Id: SmithWatermanGotoh.java,v 1.18 2005/04/03 19:38:21 ahmed Exp $
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package ccaligner;

import ccaligner.matrix.Matrix;

import java.util.ArrayList;
import java.util.logging.Logger;

import org.biojavax.bio.seq.RichSequence;

/**
 * An implementation of the Smith-Waterman algorithm with Gotoh's improvement
 * for biological local pairwise sequence alignment.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class SmithWatermanGotoh {
	/**
	 * Hidden constructor
	 */
	private SmithWatermanGotoh() {
		super();
	}

	/**
	 * Logger
	 */
	private static final Logger logger = Logger
			.getLogger(SmithWatermanGotoh.class.getName());

	/**
	 * Aligns two sequences by Smith-Waterman algorithm
	 * 
	 * @param s1
	 *            sequence #1 ({@link RichSequence})
	 * @param s2
	 *            sequence #2 ({@link RichSequence})
	 * @param pc1
	 * 			  coiled-coil registers for sequence 1 
	 * @param pc2 
	 * 			  coiled-coil registers for sequence 2 
	 * @param matrices 
	 *            coiled-coil scoring matrices ({@link Matrix})
	 * @param blosum
	 *            scoring matrix ({@link Matrix})
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @param c 
	 * @return alignment object contains the two aligned sequences, the
	 *         alignment score and alignment statistics
	 * @see RichSequence
	 * @see Matrix
	 */
	public static Alignment align(RichSequence s1, RichSequence s2, RichSequence pc1, RichSequence pc2, ArrayList<Matrix> matrices, Matrix blosum,
			float o, float e, float c_match, float c_mismatch) {
		logger.info("Started...");
		long start = System.currentTimeMillis();
		float[][] blosum_scores = blosum.getScores();
		float[][][] coil_scores = null;
		
		if (matrices != null)
		{
			coil_scores = new float[matrices.size()][][];
			
			for (int i=0; i<matrices.size(); i++)
			{
				coil_scores[i] = matrices.get(i).getScores();
		    }
		}
 
		SmithWatermanGotoh sw = new SmithWatermanGotoh();

		int m = s1.length() + 1;
		int n = s2.length() + 1;
		
		logger.info(((Integer)s1.length()).toString());
		logger.info(((Integer)pc1.length()).toString());

		
		assert s1.length() == pc1.length();
		assert s2.length() == pc2.length();

		byte[] pointers = new byte[m * n];

		// Initializes the boundaries of the traceback matrix to STOP.
		for (int i = 0, k = 0; i < m; i++, k += n) {
			pointers[k] = Directions.STOP;
		}
		for (int j = 1; j < n; j++) {
			pointers[j] = Directions.STOP;
		}

		short[] sizesOfVerticalGaps = new short[m * n];
		short[] sizesOfHorizontalGaps = new short[m * n];
		for (int i = 0, k = 0; i < m; i++, k += n) {
			for (int j = 0; j < n; j++) {
				sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
			}
		}

		Cell cell = sw.construct(s1, s2, pc1, pc2, blosum_scores, coil_scores, o, e, c_match, c_mismatch, pointers,
				sizesOfVerticalGaps, sizesOfHorizontalGaps);
		Alignment alignment = sw.traceback(s1, s2, pc1, pc2, blosum, pointers, cell,
				sizesOfVerticalGaps, sizesOfHorizontalGaps);
		alignment.setName1(s1.getName());
		alignment.setName2(s2.getName());
		alignment.setMatrix(blosum);
		alignment.setOpen(o);
		alignment.setExtend(e);
		logger.info("Finished in " + (System.currentTimeMillis() - start)
				+ " milliseconds");
		return alignment;
	}

	/**
	 * Constructs directions matrix for the traceback
	 * 
	 * @param seq1
	 *            sequence #1
	 * @param seq2
	 *            sequence #2
	 * @param pc2 
	 * @param pc1 
	 * @param blosum
	 *            scoring matrix
	 * @param coil_scores 
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @param c3 
	 * 			  coil mismatch penalty
	 * @return The cell where the traceback starts.
	 */
	private Cell construct(RichSequence seq1, RichSequence seq2, RichSequence pc1, RichSequence pc2, float[][] blosum, float[][][] coil_scores, float o,
			float e, float c_match, float c_mismatch, byte[] pointers, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		logger.info("Started...");
		long start = System.currentTimeMillis();
		
		char[] a1 = seq1.seqString().toCharArray();
		char[] a2 = seq2.seqString().toCharArray();
		int[] c1 = new int[pc1.length()];
		int[] c2 = new int[pc2.length()];

		char[] cc = pc1.seqString().toCharArray(); 
		for (int i = 0; i < pc1.length(); i++)
		{
			char c = cc[i];
			c1[i] = (c == '-') ? -1 : ((int)c - (int)'A');
		}
		
		cc = pc2.seqString().toCharArray(); 
		for (int i = 0; i < pc2.length(); i++)
		{
			char c = cc[i];
			c2[i] = (c == '-') ? -1 : ((int)c - (int)'A');
		}

		int m = seq1.length() + 1;
		int n = seq2.length() + 1;

		float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
		float[] g = new float[n]; // score if xi aligns to a gap after yi
		float h; // score if yi aligns to a gap after xi
		float[] v = new float[n]; // best score of alignment x1...xi to y1...yi
		float vDiagonal;

		g[0] = Float.NEGATIVE_INFINITY;
		h = Float.NEGATIVE_INFINITY;
		v[0] = 0;

		for (int j = 1; j < n; j++) {
			g[j] = Float.NEGATIVE_INFINITY;
			v[j] = 0;
		}

		float similarityScore, g1, g2, h1, h2;
		int r1, r2;
		char s1, s2;
		
		Cell cell = new Cell();

		for (int i = 1, k = n; i < m; i++, k += n) {
			h = Float.NEGATIVE_INFINITY;
			vDiagonal = v[0];
			for (int j = 1, l = k + 1; j < n; j++, l++) {
				
				s1 = a1[i-1];
				s2 = a2[j-1];

				if (coil_scores != null)
				{
					r1 = c1[i-1];
					r2 = c2[j-1];
					
					if (r1 >= 0 || r2 >= 0)
					{
						// use coiled-coil matrix based on query sequence
						similarityScore = coil_scores[((r1 >= 0)?r1:r2)][s1][s2]; 
						
						// at least one of the sequences is in a coil
						if (r1 == r2)
						{
							// same register: reward
							similarityScore += c_match;
						} else if (r1 >= 0 && r2 >= 0) {
							// both are coils, but in different registers: punish
							similarityScore -= c_mismatch;
						}
					}
					else
					{
						similarityScore = blosum[s1][s2];
					}
				}
				else
				{
					similarityScore = blosum[s1][s2];
				}
				
				// Fill the matrices
				f = vDiagonal + similarityScore;

				g1 = g[j] - e;
				g2 = v[j] - o;
				if (g1 > g2) {
					g[j] = g1;
					sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
				} else {
					g[j] = g2;
				}

				h1 = h - e;
				h2 = v[j - 1] - o;
				if (h1 > h2) {
					h = h1;
					sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
				} else {
					h = h2;
				}

				vDiagonal = v[j];
				v[j] = maximum(f, g[j], h, 0);

				// Determine the traceback direction
				if (v[j] == 0) {
					pointers[l] = Directions.STOP;
				} else if (v[j] == f) {
					pointers[l] = Directions.DIAGONAL;
				} else if (v[j] == g[j]) {
					pointers[l] = Directions.UP;
				} else {
					pointers[l] = Directions.LEFT;
				}

				// Set the traceback start at the current cell i, j and score
				if (v[j] > cell.getScore()) {
					cell.set(i, j, v[j]);
				}
			}
		}
		logger.info("Finished in " + (System.currentTimeMillis() - start)
				+ " milliseconds");
		return cell;
	}

	/**
	 * Returns the alignment of two sequences based on the passed array of
	 * pointers
	 * 
	 * @param s1
	 *            sequence #1
	 * @param s2
	 *            sequence #2
	 * @param pc2 
	 * @param pc1 
	 * @param m
	 *            scoring matrix
	 * @param cell
	 *            The cell where the traceback starts.
	 * @return {@link Alignment}with the two aligned sequences and alignment
	 *         score.
	 * @see Cell
	 * @see Alignment
	 */
	private Alignment traceback(RichSequence s1, RichSequence s2, RichSequence pc1, RichSequence pc2, Matrix m,
			byte[] pointers, Cell cell, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		logger.info("Started...");
		long start = System.currentTimeMillis();
		
		char[] a1 = s1.seqString().toCharArray();
		char[] a2 = s2.seqString().toCharArray();
		char[] b1 = pc1.seqString().toCharArray();
		char[] b2 = pc2.seqString().toCharArray();
		
		float[][] scores = m.getScores(); // scores at this point are only for stating similarity

		int n = s2.length() + 1;

		Alignment alignment = new Alignment();
		alignment.setScore(cell.getScore());

		int maxlen = s1.length() + s2.length(); // maximum length after the
												// aligned sequences

		char[] reversed1 = new char[maxlen]; // reversed sequence #1
		char[] reversed2 = new char[maxlen]; // reversed sequence #2
		char[] reversed3 = new char[maxlen]; // reversed markup
		char[] revcoils1 = new char[maxlen]; // reversed coils #1
		char[] revcoils2 = new char[maxlen]; // reversed coils #2

		int len1 = 0; // length of sequence #1 after alignment
		int len2 = 0; // length of sequence #2 after alignment
		int len3 = 0; // length of the markup line

		int identity = 0; // count of identitcal pairs
		int similarity = 0; // count of similar pairs
		int gaps = 0; // count of gaps

		char c1, c2;

		int i = cell.getRow(); // traceback start row
		int j = cell.getCol(); // traceback start col
		int k = i * n;

		boolean stillGoing = true; // traceback flag: true -> continue & false
								   // -> stop

		while (stillGoing) {
			switch (pointers[k + j]) {
			case Directions.UP:
				for (int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++) {
					--i;
					reversed1[len1] = a1[i];
					revcoils1[len1] = b1[i];
					len1++;
					revcoils2[len2] = Markups.GAP;
					reversed2[len2] = Alignment.GAP;
					len2++;
					reversed3[len3++] = Markups.GAP;
					k -= n;
					gaps++;
				}
				break;
				
			case Directions.DIAGONAL:
				--i;
				--j;
				c1 = a1[i];
				c2 = a2[j];
				k -= n;
				reversed1[len1] = c1;
				reversed2[len2] = c2;
				revcoils1[len1] = b1[i];
				revcoils2[len2] = b2[j];
				len1++;
				len2++;
				if (c1 == c2) {
					reversed3[len3++] = Markups.IDENTITY;
					identity++;
					similarity++;
				} else if (scores[c1][c2] > 0) {
					reversed3[len3++] = Markups.SIMILARITY;
					similarity++;
				} else {
					reversed3[len3++] = Markups.MISMATCH;
				}
				break;
				
			case Directions.LEFT:
				for (int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++) {
					reversed1[len1] = Alignment.GAP;
					revcoils1[len1] = Markups.GAP;
					len1++;
					--j;
					reversed2[len2] = a2[j];
					revcoils2[len2] = b2[j];
					len2++;
					reversed3[len3++] = Markups.GAP;
					gaps++;
				}
				break;
			case Directions.STOP:
				stillGoing = false;
			}
		}

		alignment.setSequence1(reverse(reversed1, len1));
		alignment.setStart1(i);
		alignment.setCoils1(reverse(revcoils1, len1));
		assert alignment.getSequence1().length == alignment.getCoils1().length; 
		
		alignment.setSequence2(reverse(reversed2, len2));
		alignment.setStart2(j);
		alignment.setCoils2(reverse(revcoils2, len2));
		assert alignment.getSequence2().length == alignment.getCoils2().length; 

		alignment.setMarkupLine(reverse(reversed3, len3));
		alignment.setIdentity(identity);
		alignment.setGaps(gaps);
		alignment.setSimilarity(similarity);

		logger.info("Finished in " + (System.currentTimeMillis() - start)
				+ " milliseconds");
		return alignment;
	}

	/**
	 * Returns the maximum of 4 float numbers.
	 * 
	 * @param a
	 *            float #1
	 * @param b
	 *            float #2
	 * @param c
	 *            float #3
	 * @param d
	 *            float #4
	 * @return The maximum of a, b, c and d.
	 */
	private static float maximum(float a, float b, float c, float d) {
		if (a > b) {
			if (a > c) {
				return a > d ? a : d;
			} else {
				return c > d ? c : d;
			}
		} else if (b > c) {
			return b > d ? b : d;
		} else {
			return c > d ? c : d;
		}
	}

	/**
	 * Reverses an array of chars
	 * 
	 * @param a
	 * @param len
	 * @return the input array of char reserved
	 */
	private static char[] reverse(char[] a, int len) {
		char[] b = new char[len];
		for (int i = len - 1, j = 0; i >= 0; i--, j++) {
			b[j] = a[i];
		}
		return b;
	}
}