package ccaligner.util;

import org.biojava.bio.seq.Feature.Template;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojavax.Namespace;
import org.biojavax.RankedCrossRef;
import org.biojavax.RankedDocRef;
import org.biojavax.bio.BioEntryRelationship;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.io.RichSeqIOListener;
import org.biojavax.bio.taxa.NCBITaxon;

/**
 * 
 * @author mkuhn
 *
 * Dummy implementation of the RichSeqIOListener interface from BioJava,
 * only records the sequence name.
 */
public class RichSeqIONameListener implements RichSeqIOListener {

	private String name;
	
	public void setName(String name) throws ParseException {
		this.name = name;
	}

	public String getName() {
		return name;
	}


	/**
	 * All declarations below don't do anything.
	 */
	
	public void startSequence() throws ParseException {}

	public void endSequence() throws ParseException {}

	public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length)
			throws IllegalAlphabetException {}

	public void addSequenceProperty(Object key, Object value)
			throws ParseException {}

	public void startFeature(Template templ) throws ParseException {}

	public void endFeature() throws ParseException {}

	public void addFeatureProperty(Object key, Object value)
			throws ParseException {}

	public void setAccession(String accession) throws ParseException {}

	public void setIdentifier(String identifier) throws ParseException {}

	public void setDivision(String division) throws ParseException {}

	public void setDescription(String description) throws ParseException {}

	public void setVersion(int version) throws ParseException {}

	public void setSeqVersion(String version) throws ParseException {}

	public void setComment(String comment) throws ParseException {}

	public void setRankedDocRef(RankedDocRef ref) throws ParseException {}

	public void setTaxon(NCBITaxon taxon) throws ParseException {}

	public void setNamespace(Namespace namespace) throws ParseException {}

	public void setRelationship(BioEntryRelationship relationship)
			throws ParseException {}

	public void setRankedCrossRef(RankedCrossRef crossRef)
			throws ParseException {}

	public void setURI(String uri) throws ParseException {}

	public RichFeature getCurrentFeature() throws ParseException {
		return null;
	}

	public void setCircular(boolean circular) throws ParseException {}

}
