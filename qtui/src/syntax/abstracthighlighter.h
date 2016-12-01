/*
    Copyright (C) 2016 Volker Krause <vkrause@kde.org>

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU Library General Public License as published by
    the Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
    License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __SYNTAX_ABSTRACTHIGHLIGHTER_H
#define __SYNTAX_ABSTRACTHIGHLIGHTER_H

#include <memory>
#include <boost/noncopyable.hpp>

class QString;

namespace regina {
namespace syntax {

class AbstractHighlighterPrivate;
class Definition;
class Format;
class State;
class Theme;

/**
 * Abstract base class for highlighters.
 *
 * @section abshl_intro Introduction
 *
 * The AbstractHighlighter provides an interface to highlight text.
 *
 * The SyntaxHighlighting framework already ships with one implementation,
 * namely the SyntaxHighlighter, which also derives from QSyntaxHighlighter,
 * meaning that it can be used to highlight a QTextDocument or a QML TextEdit.
 * In order to use the SyntaxHighlighter, just call setDefinition() and
 * setTheme(), and the associated documents will automatically be highlighted.
 *
 * However, if you want to use the SyntaxHighlighting framework to implement
 * your own syntax highlighter, you need to sublcass from AbstractHighlighter.
 *
 * @section abshl_impl Implementing your own Syntax Highlighter
 *
 * In order to implement your own syntax highlighter, you need to inherit from
 * AbstractHighlighter. Then, pass each text line that needs to be highlighted
 * in order to highlightLine(). Internally, highlightLine() uses the Definition
 * initially set through setDefinition() and the State of the previous text line
 * to parse and highlight the given text line. For each visual highlighting
 * change, highlightLine() will call applyFormat(). Therefore, reimplement
 * applyFormat() to get notified of the Format that is valid in the range
 * starting at the given offset with the specified length.
 *
 * The Format class itself depends on the current Theme. A theme must be
 * initially set once such that the Format%s instances can be queried for
 * concrete colors.
 *
 * Optionally, you can also reimplement setTheme() and setDefinition() to get
 * notified whenever the Definition or the Theme changes.
 *
 * @see SyntaxHighlighter
 */
class AbstractHighlighter : public boost::noncopyable
{
public:
    virtual ~AbstractHighlighter();

    /**
     * Returns the syntax definition used for highlighting.
     *
     * @see setDefinition()
     */
    Definition definition() const;

    /**
     * Sets the syntax definition used for highlighting.
     *
     * Subclasses can re-implement this method to e.g. trigger
     * re-highlighting or clear internal data structures if needed.
     */
    virtual void setDefinition(const Definition &def);

    /**
     * Returns the currently selected theme for highlighting.
     *
     * @note If no Theme was set through setTheme(), the returned Theme will be
     *       invalid, see Theme::isValid().
     */
    Theme theme() const;

    /**
     * Sets the theme used for highlighting.
     *
     * Subclasses can re-implement this method to e.g. trigger
     * re-highlighing or to do general palette color setup.
     */
    virtual void setTheme(const Theme &theme);

protected:
    AbstractHighlighter();
    AbstractHighlighter(AbstractHighlighterPrivate *dd);

    /**
     * Highlight the given line. Call this from your derived class
     * where appropriate. This will result in any number of applyFormat()
     * calls as a result.
     * @param text A string containing the text of the line to highlight.
     * @param state The highlighting state handle returned by the call
     *        to highlightLine() for the previous line. For the very first line,
     *        just pass a default constructed State().
     * @returns The state of the highlighing engine after processing the
     *        given line. This needs to passed into highlightLine() for the
     *        next line. You can store the state for efficient partial
     *        re-highlighting for example during editing.
     *
     * @see applyFormat()
     */
    State highlightLine(const QString &text, const State &state);

    /**
     * Reimplement this to apply formats to your output. The provided @p format
     * is valid for the interval [@p offset, @p offset + @p length).
     * Both @p offset and @p length count unicode characters (not machine-level
     * characters as might be used for instance in a UTF-8 encoding).
     *
     * @param offset The start column of the interval for which @p format matches
     * @param length The length of the matching text
     * @param format The Format that applies to the range [offset, offset + length)
     *
     * @note Make sure to set a valid Definition, otherwise the parameter
     *       @p format is invalid for the entire line passed to highlightLine()
     *       (cf. Format::isValid()).
     *
     * @see highlightLine()
     */
    virtual void applyFormat(int offset, int length, const Format &format) = 0;

protected:
    AbstractHighlighterPrivate *d_ptr;
};

} } // namespace regina::syntax

#endif
