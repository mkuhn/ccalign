/*
 * $Id: MatrixLoader.java,v 1.2 2005/01/25 11:54:30 ahmed Exp $
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

package ccaligner.matrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.JarURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.logging.Level;
import java.util.logging.Logger;

import ccaligner.util.Commons;


/**
 * Scoring matrices loader from a jar file or a file system.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 * @author mkuhn
 */

public class MatrixLoader {
	/**
	 * The starter character of a comment line.
	 */
	private static final char COMMENT_STARTER = '#';
	
	/**
	 * The size of the scoring matrix. It is the number of the characters in the ASCII table.
	 * It is more than the 20 amino acids just to save the processing time of the mapping. 
	 */
	private static final int SIZE = 127;
	
	/**
	 * The path to the matrices within the package.
	 */
	private static final String MATRICES_HOME = "ccaligner/matrix/matrices/";
	
	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(MatrixLoader.class.getName());

	/**
	 * List of matrix coords, to be read from header. Used to parse matrix definitions
	 */
	private final Collection<MatrixCoord> matrix_coords;
	
	/**
	 * Constructor, extract format of matrix string from header
	 * @param header
	 */
	public MatrixLoader(String header) {

		StringTokenizer tokenizer = new StringTokenizer ( header.trim(), "\t");
		
		matrix_coords = new ArrayList<MatrixCoord>(tokenizer.countTokens());

		String s = tokenizer.nextToken();
		assert(s == "# query" || s == "## query" );
		s = tokenizer.nextToken();
		assert(s == "subject");
		s = tokenizer.nextToken();
		assert(s == "scaling_factor");
		
		while (tokenizer.hasMoreTokens()) matrix_coords.add( new MatrixCoord(tokenizer.nextToken()) );
		
	}

	/**
	 * Constructor, in case there is no format line specified call the real constructor
	 * with the default format. 
	 */
	public MatrixLoader() {
		this("## query	subject	scaling_factor	" + 
				"--	-A	-B	-C	-D	-E	-F	-G	-H	-I	-K	-L	-M	-N	-P	-Q	-R	-S	-T	-V	-W	-X	-Y	-Z	-U	-*	-O	-J	" +
				"A-	AA	AB	AC	AD	AE	AF	AG	AH	AI	AK	AL	AM	AN	AP	AQ	AR	AS	AT	AV	AW	AX	AY	AZ	AU	A*	AO	AJ	" +
				"B-	BA	BB	BC	BD	BE	BF	BG	BH	BI	BK	BL	BM	BN	BP	BQ	BR	BS	BT	BV	BW	BX	BY	BZ	BU	B*	BO	BJ	" +
				"C-	CA	CB	CC	CD	CE	CF	CG	CH	CI	CK	CL	CM	CN	CP	CQ	CR	CS	CT	CV	CW	CX	CY	CZ	CU	C*	CO	CJ	" +
				"D-	DA	DB	DC	DD	DE	DF	DG	DH	DI	DK	DL	DM	DN	DP	DQ	DR	DS	DT	DV	DW	DX	DY	DZ	DU	D*	DO	DJ	" +
				"E-	EA	EB	EC	ED	EE	EF	EG	EH	EI	EK	EL	EM	EN	EP	EQ	ER	ES	ET	EV	EW	EX	EY	EZ	EU	E*	EO	EJ	" +
				"F-	FA	FB	FC	FD	FE	FF	FG	FH	FI	FK	FL	FM	FN	FP	FQ	FR	FS	FT	FV	FW	FX	FY	FZ	FU	F*	FO	FJ	" +
				"G-	GA	GB	GC	GD	GE	GF	GG	GH	GI	GK	GL	GM	GN	GP	GQ	GR	GS	GT	GV	GW	GX	GY	GZ	GU	G*	GO	GJ	" +
				"H-	HA	HB	HC	HD	HE	HF	HG	HH	HI	HK	HL	HM	HN	HP	HQ	HR	HS	HT	HV	HW	HX	HY	HZ	HU	H*	HO	HJ	" +
				"I-	IA	IB	IC	ID	IE	IF	IG	IH	II	IK	IL	IM	IN	IP	IQ	IR	IS	IT	IV	IW	IX	IY	IZ	IU	I*	IO	IJ	" +
				"K-	KA	KB	KC	KD	KE	KF	KG	KH	KI	KK	KL	KM	KN	KP	KQ	KR	KS	KT	KV	KW	KX	KY	KZ	KU	K*	KO	KJ	" +
				"L-	LA	LB	LC	LD	LE	LF	LG	LH	LI	LK	LL	LM	LN	LP	LQ	LR	LS	LT	LV	LW	LX	LY	LZ	LU	L*	LO	LJ	" +
				"M-	MA	MB	MC	MD	ME	MF	MG	MH	MI	MK	ML	MM	MN	MP	MQ	MR	MS	MT	MV	MW	MX	MY	MZ	MU	M*	MO	MJ	" +
				"N-	NA	NB	NC	ND	NE	NF	NG	NH	NI	NK	NL	NM	NN	NP	NQ	NR	NS	NT	NV	NW	NX	NY	NZ	NU	N*	NO	NJ	" +
				"P-	PA	PB	PC	PD	PE	PF	PG	PH	PI	PK	PL	PM	PN	PP	PQ	PR	PS	PT	PV	PW	PX	PY	PZ	PU	P*	PO	PJ	" +
				"Q-	QA	QB	QC	QD	QE	QF	QG	QH	QI	QK	QL	QM	QN	QP	QQ	QR	QS	QT	QV	QW	QX	QY	QZ	QU	Q*	QO	QJ	" +
				"R-	RA	RB	RC	RD	RE	RF	RG	RH	RI	RK	RL	RM	RN	RP	RQ	RR	RS	RT	RV	RW	RX	RY	RZ	RU	R*	RO	RJ	" +
				"S-	SA	SB	SC	SD	SE	SF	SG	SH	SI	SK	SL	SM	SN	SP	SQ	SR	SS	ST	SV	SW	SX	SY	SZ	SU	S*	SO	SJ	" +
				"T-	TA	TB	TC	TD	TE	TF	TG	TH	TI	TK	TL	TM	TN	TP	TQ	TR	TS	TT	TV	TW	TX	TY	TZ	TU	T*	TO	TJ	" +
				"V-	VA	VB	VC	VD	VE	VF	VG	VH	VI	VK	VL	VM	VN	VP	VQ	VR	VS	VT	VV	VW	VX	VY	VZ	VU	V*	VO	VJ	" +
				"W-	WA	WB	WC	WD	WE	WF	WG	WH	WI	WK	WL	WM	WN	WP	WQ	WR	WS	WT	WV	WW	WX	WY	WZ	WU	W*	WO	WJ	" +
				"X-	XA	XB	XC	XD	XE	XF	XG	XH	XI	XK	XL	XM	XN	XP	XQ	XR	XS	XT	XV	XW	XX	XY	XZ	XU	X*	XO	XJ	" +
				"Y-	YA	YB	YC	YD	YE	YF	YG	YH	YI	YK	YL	YM	YN	YP	YQ	YR	YS	YT	YV	YW	YX	YY	YZ	YU	Y*	YO	YJ	" +
				"Z-	ZA	ZB	ZC	ZD	ZE	ZF	ZG	ZH	ZI	ZK	ZL	ZM	ZN	ZP	ZQ	ZR	ZS	ZT	ZV	ZW	ZX	ZY	ZZ	ZU	Z*	ZO	ZJ	" +
				"U-	UA	UB	UC	UD	UE	UF	UG	UH	UI	UK	UL	UM	UN	UP	UQ	UR	US	UT	UV	UW	UX	UY	UZ	UU	U*	UO	UJ	" +
				"*-	*A	*B	*C	*D	*E	*F	*G	*H	*I	*K	*L	*M	*N	*P	*Q	*R	*S	*T	*V	*W	*X	*Y	*Z	*U	**	*O	*J	" +
				"O-	OA	OB	OC	OD	OE	OF	OG	OH	OI	OK	OL	OM	ON	OP	OQ	OR	OS	OT	OV	OW	OX	OY	OZ	OU	O*	OO	OJ	" +
				"J-	JA	JB	JC	JD	JE	JF	JG	JH	JI	JK	JL	JM	JN	JP	JQ	JR	JS	JT	JV	JW	JX	JY	JZ	JU	J*	JO	JJ"
			);
	}

	/**
	 * Read matrix from line, using the header definition read by the constructor
	 * @param line
	 * @return
	 */
	public Matrix loadFromLine(String line)
	{
		float[][] scores = new float[SIZE][SIZE];
		
		if (line.startsWith("#")) line = line.substring(1);
		
		StringTokenizer tokenizer = new StringTokenizer ( line.trim(), "\t");

		if (tokenizer.countTokens() < matrix_coords.size()) return null;
		
		String query = tokenizer.nextToken();
		String subject = tokenizer.nextToken();
		Float scaling_factor = Float.valueOf(tokenizer.nextToken());
		
		Iterator<MatrixCoord> it = matrix_coords.iterator();
		
		while (tokenizer.hasMoreTokens())
		{
			MatrixCoord mc = it.next();
			String t = tokenizer.nextToken();
			Float v = (float) -32768;
			if (!t.contentEquals("X")) v = Float.valueOf(t) / scaling_factor;
			scores[mc.getI()][mc.getJ()] = v; 
		}
		
		return new Matrix(query, subject, scores);
	}
	
	/**
	 * Loads scoring matrices from Jar file or file system.
	 * @param matrix to load
	 * @return loaded matrices as map
	 * @throws MatrixLoaderException
	 * @see Matrix
	 */
	public static Map<String,Matrix> loadMatrices (String matrix) throws MatrixLoaderException {
	    logger.fine("Trying to load scoring matrices... " + matrix );

		InputStream is = null;
		
		if (new StringTokenizer(matrix, Commons.getFileSeparator()).countTokens() == 1) {
			// Matrix does not include the path
			// Load the matrix from matrices.jar
			is = MatrixLoader.class.getClassLoader().getResourceAsStream(MATRICES_HOME + matrix);
			
			if (is == null)
			{
		        String message = "Failed opening input stream! " + MATRICES_HOME + matrix;
		        logger.log(Level.SEVERE, message);
		        throw new MatrixLoaderException (message);
			}
		} else {
			// Matrix includes the path information
			// Load the matrix from the file system
			try {
			    is = new FileInputStream(matrix);
		    } catch (Exception e) {
		        String message = "Failed opening input stream: " + e.getMessage();
		        logger.log(Level.SEVERE, message, e);
		        throw new MatrixLoaderException (message);
		    }
		}

		return loadMatrices(matrix, is);
	}	
	
	/**
	 * Loads scoring matrix from {@link InputStream}
	 * @param nis named input stream
	 * @return loaded matrix
	 * @throws MatrixLoaderException
	 * @see Matrix
	 * @see NamedInputStream
	 */
	public static Map<String,Matrix> loadMatrices (String matrix, InputStream is) throws MatrixLoaderException {
	    logger.fine("Loading scoring matrices... " + matrix );
			
		BufferedReader reader = new BufferedReader(new InputStreamReader(is));
		
		Map<String,Matrix> matrices = new HashMap<String,Matrix>();
		
		try {
			String line = reader.readLine();

			MatrixLoader ml = new MatrixLoader(line);
			
			while ((line = reader.readLine()) != null)
			{
				Matrix m = ml.loadFromLine(line);
				matrices.put(m.getId(), m);
			}

		} catch (IOException e) {
	        String message = "Failed reading from input stream: " + e.getMessage();
	        logger.log(Level.SEVERE, message, e);
	        throw new MatrixLoaderException (message);
		}

		return matrices;
	}
	
	/**
	 * Loads scoring matrix from Jar file or file system.
	 * @param matrix to load
	 * @return loaded matrix
	 * @throws MatrixLoaderException
	 * @see Matrix
	 */
	public static Matrix load (String matrix) throws MatrixLoaderException {
		InputStream is = null;
		
		if (new StringTokenizer(matrix, Commons.getFileSeparator()).countTokens() == 1) {
			// Matrix does not include the path
			// Load the matrix from matrices.jar
			is = MatrixLoader.class.getClassLoader().getResourceAsStream(MATRICES_HOME + matrix);
		} else {
			// Matrix includes the path information
			// Load the matrix from the file system
			try {
			    is = new FileInputStream(matrix);
		    } catch (Exception e) {
		        String message = "Failed opening input stream: " + e.getMessage();
		        logger.log(Level.SEVERE, message, e);
		        throw new MatrixLoaderException (message);
		    }
		}

		return load(matrix, is);
	}

	/**
	 * Loads scoring matrix from {@link InputStream}
	 * @param nis named input stream
	 * @return loaded matrix
	 * @throws MatrixLoaderException
	 * @see Matrix
	 * @see NamedInputStream
	 */
	public static Matrix load (String matrix, InputStream is) throws MatrixLoaderException {
	    logger.fine("Loading scoring matrix... " + matrix );
	    char[] acids = new char[SIZE];
			
		// Initialize the acids array to null values (ascii = 0)
		for (int i = 0; i < SIZE; i++) {
			acids[i] = 0;
		}
			
		float[][] scores = new float[SIZE][SIZE];

		BufferedReader reader = new BufferedReader(new InputStreamReader(is));
		
		String line;
			
		try {
			// Skip the comment lines
			while ((line = reader.readLine()) != null && line.trim().charAt(0) == COMMENT_STARTER);
	    } catch (Exception e) {
	        String message = "Failed reading from input stream: " + e.getMessage();
	        logger.log(Level.SEVERE, message, e);
	        throw new MatrixLoaderException (message);
	    }
	
		// Read the headers line (the letters of the acids)
		StringTokenizer tokenizer;
		tokenizer = new StringTokenizer ( line.trim( ) );
		for (int j = 0; tokenizer.hasMoreTokens(); j++) {
			acids[j] = tokenizer.nextToken().charAt(0);
		}

		try {
			int current_acid = 0; 
			
			// Read the scores
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				tokenizer = new StringTokenizer ( line );

				// check if first char of matrix is the amino acids
				char acid = line.charAt(0);
				
				if (Character.isLetter(acid) || acid == '*')
				{
					tokenizer.nextToken();
				}
				else
				{
					// if not the case, use same order as specified in first row
					acid = acids[current_acid++];
				}
				
				for (int i = 0; i < SIZE; i++) {
					if (acids[i] != 0 && tokenizer.hasMoreTokens()) {
						scores[acids[i]][acid] = scores[acid][acids[i]] = Float.parseFloat(tokenizer.nextToken()); 
					}
				}
			}
	    } catch (Exception e) {
	        String message = "Failed reading from input stream: " + e.getMessage();
	        logger.log(Level.SEVERE, message, e);
	        throw new MatrixLoaderException (message);
	    }
	    logger.fine("Finished loading scoring matrix");
		return new Matrix(matrix, scores);
	}

	/**
	 * Returns a list of the scoring matrices in the matrices home directory
	 * @param sort flag to sort the list or not
	 * @return sorted array of scoring matrices
	 * @throws MatrixLoaderException
	 */
	public static Collection<String> list (boolean sort ) throws MatrixLoaderException {
		logger.fine("Loading list of scoring matrices...");
	    ArrayList<String> matrices = new ArrayList<String>();
		URL url = MatrixLoader.class.getClassLoader().getResource(MATRICES_HOME);
		if (url.getFile().toString().indexOf("!") != -1) {
			// Load from Jar
		    JarURLConnection connection = null;
		    JarFile jar = null;
		    try {
		        connection = (JarURLConnection) url.openConnection();
		        jar = connection.getJarFile();
		    } catch (Exception e) {
		        String message = "Failed opening a connection to jar: " + e.getMessage();
		        logger.log(Level.SEVERE, message, e);
		        throw new MatrixLoaderException (message);
		    }
			Enumeration<?> entries = jar.entries();
			JarEntry entry;
			String entryName;
			int length = MATRICES_HOME.length();
			while (entries.hasMoreElements()) {
				entry = (JarEntry) entries.nextElement();
				if (!entry.isDirectory()) {
					entryName = entry.getName();
					if (entryName.startsWith(MATRICES_HOME)) {
						matrices.add(entryName.substring(length));
					}
				}
			}
		} else {
			// Load from file system
			String home = url.getFile( );
			File dir = new File (home);
			String files[] = dir.list( );
			File file;
			for (int i = 0, n = files.length; i < n; i++) {
				file = new File (home + files[i]);
				if (file.isFile() && file.canRead()) {
					matrices.add(file.getName());
				}
			}
		}
		if (sort) {
		    Collections.sort(matrices, new MatricesComparator());
		}
		logger.fine("Finished loading list of scoring matrices");
		return matrices;
	}

	/**
	 * Returns a list of the scoring matrices in the matrices home directory
	 * @return sorted array of scoring matrices
	 * @throws MatrixLoaderException
	 */
	public static Collection<String> list ( ) throws MatrixLoaderException {
	    return list(false);
	}
}