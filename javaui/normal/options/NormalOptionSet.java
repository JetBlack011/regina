
/**************************************************************************
 *                                                                        *
 *  Regina - A normal surface theory calculator                           *
 *  Java user interface                                                   *
 *                                                                        *
 *  Copyright (c) 1999-2001, Ben Burton                                   *
 *  For further details contact Ben Burton (benb@acm.org).                *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,        *
 *  MA 02111-1307, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

package normal.options;

import java.io.*;
import java.util.Vector;
import btools.utilities.OptionSet;

/**
 * Provides an option set with a cache of frequently used options for
 * faster access.  These are accessed through direct variable lookup
 * as opposed to the usual dictionary search with string comparisons.
 * <p>
 * These frequently used options should <i>only</i> be accessed through
 * the appropriate <tt>get</tt><i>OptionName</i><tt>()</tt>
 * and <tt>set</tt><i>OptionName</i><tt>()</tt> methods, and never through
 * the standard methods <tt>getStringOption()</tt>,
 * <tt>setBooleanOption()</tt>, etc.
 * <p>
 * Options that are not explicitly cached by <tt>NormalOptionSet</tt> should,
 * however, still be accessed through the standard <tt>OptionSet</tt> access
 * methods.
 * <p>
 * Cached options will always provide values.  If they do not appear in the
 * option set, they will be given default values upon construction of the
 * option set.
 */
public class NormalOptionSet extends OptionSet {
    /**
     * The comment to place at the beginning of the option file.
     */
    public static final String normalComment = "User options for " +
        normal.Application.program;
    
    /**
     * Full option name for a particular cached option.
     */
    public static final String optionAutoDock = "AutoDock";
    
    /**
     * Full option name for a particular cached option.
     */
    public static final String optionDisplayIcon = "DisplayIcon";

    /**
     * Full option name for a particular cached option.
     */
    public static final String optionJPythonLibCount = "JPythonLibCount";

    /**
     * Option name stub for a particular cached option.
     */
    public static final String optionJPythonLib = "JPythonLib";

    /**
     * Option name stub for a particular cached option.
     */
    public static final String optionJPythonLibUse = "JPythonLibUse";

    /**
     * The default for a particular cached option.
     */
    public static final boolean defaultAutoDock = true;
    
    /**
     * The default for a particular cached option.
     */
    public static final boolean defaultDisplayIcon = true;
    
    /**
     * A particular cached option.
     */
    private boolean autoDock;
    
    /**
     * A particular cached option.
     */
    private boolean displayIcon;

	/**
	 * A particular cached option.
	 */
	private Vector jpythonLibraries;

    /**
     * Creates a new option set not associated with any file.
     */
    public NormalOptionSet() {
        super();
    }
    
    /**
     * Creates a new option set based on the given file.
     * All options stored in the file will be loaded.  Any errors
     * that occur when reading from file will be ignored.
     *
     * @param optionFile the file to/from which this option set is to
     * be written/read.
     * @see #NormalOptionSet(File, boolean)
     */
    public NormalOptionSet(File optionFile) {
        super(optionFile, normalComment);
    }

    /**
     * Creates a new option set based on the given file.
     * All options stored in the file will be loaded.
     *
     * @param optionFile the file to/from which this option set is to
     * be written/read.
     * @param forceLoad if set to <tt>false</tt>, we will ignore
     * any errors in reading from file.  If set to <tt>true</tt>, an
     * exception will be thrown if an error occurs.
     * @throws IOException thrown if an error occurs in reading from file.
     * @see #NormalOptionSet(File)
     */
    public NormalOptionSet(File optionFile, boolean forceLoad)
            throws IOException {
        super(optionFile, normalComment, forceLoad);
    }
    
    /**
     * Attempt to read this option set from file.
     *
     * @throws IOException thrown if an error occurs in reading from file.
     */
    protected void readFromFile() throws IOException {
        IOException caught = null;
        try {
            super.readFromFile();
        } catch (IOException e) {
            caught = e;
        }
        autoDock = getBooleanOption(optionAutoDock, defaultAutoDock);
        displayIcon = getBooleanOption(optionDisplayIcon, defaultDisplayIcon);

		jpythonLibraries = new Vector();
		int nLibs = getIntOption(optionJPythonLibCount, 0);
		for (int i = 0; i < nLibs; i++)
			jpythonLibraries.addElement(new JPythonLibrary(
				getStringOption(optionJPythonLib + String.valueOf(i), ""),
				getBooleanOption(optionJPythonLibUse + String.valueOf(i),
				false)));

        if (caught != null)
            throw caught;
    }

    /**
     * Attempt to write this option set to file.
     *
     * @param forceWrite if set to <tt>false</tt>, we will ignore
     * any errors in writing to file.  If set to <tt>true</tt>, an
     * exception will be thrown if an error occurs.
     * @throws IOException thrown if an error occurs in writing to file.
     */
    public void writeToFile(boolean forceWrite) throws IOException {
        setBooleanOption(optionAutoDock, autoDock);
        setBooleanOption(optionDisplayIcon, displayIcon);

		int nLibs = jpythonLibraries.size();
		// Clear out all unnecessary library options.
		int oldNLibs = getIntOption(optionJPythonLibCount, 0);
		for (int i = nLibs; i < oldNLibs; i++) {
			removeOption(optionJPythonLib + String.valueOf(i));
			removeOption(optionJPythonLibUse + String.valueOf(i));
		}

		setIntOption(optionJPythonLibCount, nLibs);
		JPythonLibrary lib;
		for (int i = 0; i < nLibs; i++) {
			lib = (JPythonLibrary)jpythonLibraries.elementAt(i);
			setStringOption(optionJPythonLib + String.valueOf(i),
				lib.getLibraryPath());
			setBooleanOption(optionJPythonLibUse + String.valueOf(i),
				lib.shouldUseLibrary());
		}

        super.writeToFile(forceWrite);
    }

    /**
     * Get a particular cached option.
     *
     * @return the current value of the option.
     */
    public boolean getAutoDock() {
        return autoDock;
    }
    
    /**
     * Set a particular cached option.
     *
     * @param value the new value for the option.
     */
    public void setAutoDock(boolean value) {
        autoDock = value;
    }
        
    /**
     * Get a particular cached option.
     *
     * @return the current value of the option.
     */
    public boolean getDisplayIcon() {
        return displayIcon;
    }
    
    /**
     * Set a particular cached option.
     *
     * @param value the new value for the option.
     */
    public void setDisplayIcon(boolean value) {
        displayIcon = value;
    }

	/**
	 * Returns a vector of <tt>JPythonLibrary</tt> objects representing
	 * the set of available JPython libraries.
	 *
	 * @return the set of available JPython libraries.
	 * @see normal.options.NormalOptionSet.JPythonLibrary
	 */
	public Vector getJPythonLibraries() {
		return jpythonLibraries;
	}

	/**
	 * Sets the list of available JPython libraries.
	 *
	 * @param value the new list of JPython libraries; this must be a
	 * vector of <tt>JPythonLibrary</tt> objects.
	 * @see normal.options.NormalOptionSet.JPythonLibrary
	 */
	public void setJPythonLibraries(Vector value) {
		jpythonLibraries = value;
	}

	/**
	 * Stores the details of a JPython library.
	 */
	public static class JPythonLibrary {
		/**
		 * The path to the JPython library file.
		 */
		private String libraryPath;
		/**
		 * Specifies whether or not the library should be used.
		 */
		private boolean useLibrary;

		/**
		 * Creates a new JPython library specifier.
		 *
		 * @param libraryPath the path to the JPython library file.
		 * @param useLibrary <tt>true</tt> if and only if the library
		 * should be used.
		 */
		public JPythonLibrary(String libraryPath, boolean useLibrary) {
			this.libraryPath = libraryPath;
			this.useLibrary = useLibrary;
		}

		/**
		 * Returns the path to this particular JPython library file.
		 *
		 * @return the path to this particular JPython library file.
		 */
		public String getLibraryPath() {
			return libraryPath;
		}

		/**
		 * Returns whether or not this JPython library should be used.
		 *
		 * @return <tt>true</tt> if and only if this JPython library
		 * should be used.
		 */
		public boolean shouldUseLibrary() {
			return useLibrary;
		}
	}
}

