
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

package normal.console;

import java.io.*;
import java.util.*;
import normal.Shell;
import normal.options.*;
import org.python.core.*;
import org.python.modules.*;
import org.python.util.*;

/**
 * Provides various utilities for compiling and running JPython code.
 */
public class JPythonUtils {
	/**
	 * The commands to be run whenever a new JPython session is started.
	 * These commands will be run (in the order in which they appear in
	 * the array) before anything else is done.
	 * Each command should be a single element in the array.
	 */
	public static final String[] startup = {
		"import java", "from java.lang import *",
		"import normal", "import btools"
	};

	/**
	 * Initialises the given JPython interpreter.
	 * The startup commands in <tt>JPythonUtils.startup</tt> are run
	 * and the variable <tt>engine</tt> is set to the given calculation
	 * engine.
	 *
	 * @param interpreter the JPython interpreter to set up; this should
	 * be of class <tt>org.python.util.PythonInterpreter</tt> but is
	 * passed as an <tt>Object</tt> to ease stability in the cases where
	 * the python classes are not available.
	 * @param shell the shell representing the entire program.
	 * @return a text string to inform the user of the initialisation
	 * that has been done; this may contain multiple lines and will end
	 * in a final newline.
	 */
	public static String setupInterpreter(Object interpreter, Shell shell) {
		PythonInterpreter realInterpreter = (PythonInterpreter)interpreter;
		PyObject code;
		StringBuffer error = new StringBuffer();
		String message;

		message = "Running startup commands.\n";
		for (int i = 0; i < startup.length; i++) {
			code = compileCode(startup[i], error);
			if (code == null) {
				message = message + "Error compiling: " + startup[i] +
					'\n' + error + '\n';
				error = new StringBuffer();
			} else {
				if (! runCode(code, realInterpreter, error)) {
					message = message + "Error running: " + startup[i] +
						'\n' + error + '\n';
					error = new StringBuffer();
				}
			}
		}

		// Set the engine.
		realInterpreter.set("engine", shell.getEngine());
		message = message + "The calculation engine " +
			"(type normal.engine.Engine) is in the variable [engine].\n";

		// Load libraries.
		File libFile;
		NormalOptionSet.JPythonLibrary lib;
		Enumeration e = shell.getOptions().getJPythonLibraries().elements();
		while (e.hasMoreElements()) {
			lib = (NormalOptionSet.JPythonLibrary)e.nextElement();
			if (lib.shouldUseLibrary()) {
				libFile = new File(lib.getLibraryPath());
				message = message + "Loading: " + libFile.getName() + '\n';

				code = compileFile(libFile, error);
				if (code == null) {
					message = message + error + '\n';
					error = new StringBuffer();
				} else {
					if (! runCode(code, realInterpreter, error)) {
						message = message + error + '\n';
						error = new StringBuffer();
					}
				}
			}
		}

		return message;
	}

	/**
	 * Compiles the given block of JPython code.  The block of code
	 * should be presented as a number of lines of code separated by
	 * newlines.
	 * <p>
	 * If a compile error occurs, the corresponding error message will
	 * be appended to the given string buffer.
	 *
	 * @param code the block of code to compile.
	 * @param error a string buffer to which any compile errors will be
	 * appended.  Note that compile errors may contain many lines of
	 * information and will always end in a final newline.
	 * @return the compiled code on success, or <tt>null</tt> on
	 * failure.
	 */
	public static PyObject compileCode(String code, StringBuffer error) {
		PyObject compiled = null;
		try {
			compiled = codeop.compile_command(
				code, "<script>", "exec");
		} catch (RuntimeException exc) {
			if (exc instanceof PyException) {
				if (Py.matchException((PyException)exc, Py.SyntaxError)) {
					compiled = null;
					error.append(exc.toString());
				} else
					throw exc;
			} else
				throw exc;
		}
		if (compiled == Py.None) {
			compiled = null;
			error.append(
				"The script was incomplete; JPython expects more code.");
		}

		// Append a final newline to the error message if necessary.
		if (compiled == null) {
			int len = error.length();
			if (len == 0)
				error.append('\n');
			else if (error.charAt(len - 1) != '\n')
				error.append('\n');
		}

		return compiled;
	}

	/**
	 * Compiles the given JPython file.
	 * If a compile or file I/O error occurs, the corresponding error
	 * message will be appended to the given string buffer.
	 *
	 * @param file the JPython file to compile.
	 * @param error a string buffer to which any compile errors will be
	 * appended.  Note that compile errors may contain many lines of
	 * information and will always end in a final newline.
	 * @return the compiled code on success, or <tt>null</tt> on
	 * failure.
	 */
	public static PyObject compileFile(File file, StringBuffer error) {
		FileInputStream in;
		try {
			in = new FileInputStream(file);
		} catch (Exception e) {
			error.append("Could not open file [" +
				file.getAbsolutePath() + "] for reading.\n");
			error.append(e.toString());
			if (error.charAt(error.length() - 1) != '\n')
				error.append('\n');
			return null;
		}

		PyObject compiled = null;
		try {
			compiled = Py.compile(in, file.getName(), "exec");
		} catch (RuntimeException exc) {
			if (exc instanceof PyException) {
				if (Py.matchException((PyException)exc, Py.SyntaxError)) {
					compiled = null;
					error.append(exc.toString());
				} else
					throw exc;
			} else
				throw exc;
		}
		if (compiled == Py.None) {
			compiled = null;
			error.append(
				"The file was incomplete; JPython expects more code.");
		}

		// Append a final newline to the error message if necessary.
		if (compiled == null) {
			int len = error.length();
			if (len == 0)
				error.append('\n');
			else if (error.charAt(len - 1) != '\n')
				error.append('\n');
		}

		return compiled;
	}

	/**
	 * Runs the given pre-compiled JPython code.
	 * If a runtime error occurs, the corresponding error
	 * message will be appended to the given string buffer.
	 *
	 * @param compiledCode the pre-compiled JPython code to run.
	 * @param interpreter the JPython interpreter in which to run the code.
	 * @param error a string buffer to which any runtime errors will be
	 * appended.  Note that runtime errors may contain many lines of
	 * information and will always end in a final newline.
	 * @return <tt>true</tt> on success or <tt>false</tt> on failure.
	 */
	public static boolean runCode(PyObject compiledCode,
			PythonInterpreter interpreter, StringBuffer error) {
		try {
			interpreter.exec(compiledCode);
			return true;
		} catch (RuntimeException exc) {
			if (exc instanceof PyException) {
				// Runtime error!
				error.append(exc.toString());
				int len = error.length();
				if (len == 0)
					error.append('\n');
				else if (error.charAt(len - 1) != '\n')
					error.append('\n');
				return false;
			} else
				throw exc;
		}
	}
}

