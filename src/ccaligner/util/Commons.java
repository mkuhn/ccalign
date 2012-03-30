/*
 * $Id: Commons.java,v 1.27 2005/04/18 05:37:52 ahmed Exp $
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

package ccaligner.util;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Global constants/variables/settings
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public abstract class Commons {
	
	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(Commons.class
			.getName());

	/**
	 * Default home directory
	 */
	private static final String DEFAULT_USER_DIRECTORY = ".";

	/**
	 * Default file separator
	 */
	private static final String DEFAULT_FILE_SEPARATOR = "/";

	/**
	 * Default line separator
	 */
	private static final String DEFAULT_LINE_SEPARATOR = "\r\n";

	/**
	 * User home directory
	 */
	private static String userDirectory = DEFAULT_USER_DIRECTORY;
	static {
		try {
			userDirectory = System.getProperty("user.dir");
		} catch (Exception e) {
			logger.log(Level.WARNING, "Failed getting user current directory: "
					+ e.getMessage(), e);
		}
	}

	/**
	 * Line separator
	 */
	private static String fileSeparator = DEFAULT_FILE_SEPARATOR;
	static {
		try {
			fileSeparator = System.getProperty("file.separator");
		} catch (Exception e) {
			logger.log(Level.WARNING, "Failed getting system file separator: "
					+ e.getMessage(), e);
		}
	}

	/**
	 * Line separator
	 */
	private static String lineSeparator = DEFAULT_LINE_SEPARATOR;

	static {
		try {
			lineSeparator = System.getProperty("line.separator");
		} catch (Exception e) {
			logger.log(Level.WARNING, "Failed getting system line separator: "
					+ e.getMessage(), e);
		}
	}

	/**
	 * JNLP enabled
	 */
	private static boolean jnlp = false;
	static {
		try {
			jnlp = "true".equalsIgnoreCase(System.getProperty("jnlp.enabled"));
		} catch (Exception e) {
			logger.log(Level.WARNING, "Failed getting jnlp enabled property: "
					+ e.getMessage(), e);
		}
		setJnlp(jnlp);
	}

	/**
	 * Returns system file separator.
	 * 
	 * @return file separator
	 */
	public static String getFileSeparator() {
		return fileSeparator;
	}

	/**
	 * Returns system line separator.
	 * 
	 * @return line separator
	 */
	public static String getLineSeparator() {
		return lineSeparator;
	}

	/**
	 * Returns user's current directory.
	 * 
	 * @return user's current directory
	 */
	public static String getUserDirectory() {
		return userDirectory;
	}

	/**
	 * Sets the JNLP flag to true of false
	 * 
	 * @param jnlp
	 *            true or false
	 */
	public static void setJnlp(boolean jnlp) {
		Commons.jnlp = jnlp;
	}

	/**
	 * Returns true if jnlp is enabled
	 * 
	 * @return true if jnlp is enabled, otherwise returns false
	 */
	public static boolean isJnlp() {
		return jnlp;
	}
}