
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2014, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

/*! \file treewidth/ntreedecomposition.h
 *  \brief Deals with tree decompositions of facet pairing graphs.
 */

#ifndef __NTREEDECOMPOSITION_H
#ifndef __DOXYGEN
#define __NTREEDECOMPOSITION_H
#endif

#include "regina-core.h"
#include "shareableobject.h"
#include "census/ngenericfacetpairing.h"
#include "generic/ngenerictriangulation.h"

namespace regina {

/**
 * \addtogroup treewidth Treewidth
 * Treewidth and tree decompositions.
 * @{
 */

class NTreeBag;

/**
 * TODO: Document this file.
 */

enum TreeDecompositionAlg {
    TD_UPPER = 0x0001,
    TD_UPPER_GREEDY_FILL_IN = 0x0001,
    /**
     * Not yet implemented: will fall back to greedy for now.
     */
    TD_EXACT = 0x0010
};

enum BagComparison {
    BAG_EQUAL = 0,
    BAG_SUBSET = -1,
    BAG_SUPERSET = 1,
    BAG_UNRELATED = 2
};

enum NiceType {
    NICE_INTRODUCE = 1,
    NICE_FORGET = 2,
    NICE_JOIN = 3
};

class REGINA_API NTreeBag : public ShareableObject {
    private:
        int size_;
            /**< The number of elements in this bag. */
        int* elements_;
            /**< The elements of this bag, sorted in ascending order. */
        NTreeBag* parent_;
        NTreeBag* sibling_;
        NTreeBag* children_;
        int type_;
            /**< Zero if nothing special; otherwise a non-zero type
                 constant specific to the application. */
        int subtype_;
        int index_;
            /**< Undefined until set by NTreeDecomposition.
                 Otherwise ordered from leaves to root. */

    public:
        /**
         * Note: only the list of elements will be cloned.
         * The bag will not be inserted into the tree (so parent_, sibling_
         * and children_ will all be null), and type and subtype will not
         * be set.
         */
        NTreeBag(const NTreeBag& cloneMe);
        ~NTreeBag();

        int size() const;
        int element(int which) const;
        bool contains(int element) const;
        int index() const;

        int type() const;
        int subtype() const;

        BagComparison compare(const NTreeBag& rhs) const;

        const NTreeBag* next() const;
        const NTreeBag* nextPrefix() const;

        const NTreeBag* parent() const;
        const NTreeBag* children() const;
        const NTreeBag* sibling() const;
        bool isLeaf() const;

        void writeTextShort(std::ostream& out) const;

    private:
        NTreeBag(int size);

        /**
         * Inserts as the first child.
         */
        void insertChild(NTreeBag* child);

        /**
         * Only swaps the lists of elements.
         */
        void swapContents(NTreeBag& other);

    friend class NTreeDecomposition;
};

class REGINA_API NTreeDecomposition : public ShareableObject {
    protected:
        /**
         * Note: loops are ignored.
         */
        struct Graph {
            int order_;
            bool** adj_;

            /**
             * Fills with false.
             */
            Graph(int order);
            ~Graph();

            void dump(std::ostream& out) const;
        };

    private:
        int width_;
            /**< The width of the tree decomposition; that is, one less
                 than the maximum bag size. */
        int size_;
        NTreeBag* root_;

    public:
        /**
         * \ifacespython The first argument must be of type NTriangulation
         * or Dim2Triangulation.
         */
        template <int dim>
        NTreeDecomposition(
            const NGenericTriangulation<dim>& triangulation,
            TreeDecompositionAlg alg = TD_UPPER);

        /**
         * \ifacespython The first argument must be of type NFacePairing
         * or Dim2EdgePairing.
         */
        template <int dim>
        NTreeDecomposition(
            const NGenericFacetPairing<dim>& pairing,
            TreeDecompositionAlg alg = TD_UPPER);

        /**
         * Note: if the matrix is asymmetric (a digraph), then the
         * undirected graph will be used.  Loops are ignored.
         */
        template <typename T>
        NTreeDecomposition(unsigned order, T const** const graph,
            TreeDecompositionAlg alg = TD_UPPER);

        ~NTreeDecomposition();

        int width() const;
        int size() const;

        const NTreeBag* root() const;
        const NTreeBag* first() const;
        const NTreeBag* firstPrefix() const;

        /**
         * Merge adjacent bags where one is a subset of another.
         */
        bool compress();
        void makeNice();

        void writeTextShort(std::ostream& out) const;
        void writeTextLong(std::ostream& out) const;
#if 0
        /**
         * Verifies that (i) all bag elements are in range;
         * (ii) all elements appear in some bag; and
         * (iii) the path condition holds.
         * Does not verify the edge condition (since we do not have
         * access to the edges of the underlying graph).
         */
        bool verify(int order) const;

        template <int dim>
        bool verify(const NGenericTriangulation<dim>& triangulation) const;
#endif
    private:
        /**
         * Note: graph may be modified during this routine.
         */
        void construct(Graph& graph, TreeDecompositionAlg alg);
        void greedyFillIn(Graph& graph);
        void reindex();
};

/*@}*/

// Inline functions for NTreeBag

inline NTreeBag::NTreeBag(int size) :
        size_(size),
        elements_(new int[size_]),
        parent_(0),
        sibling_(0),
        children_(0),
        type_(0),
        subtype_(0) {
}

inline NTreeBag::NTreeBag(const NTreeBag& cloneMe) :
        size_(cloneMe.size_),
        elements_(new int[cloneMe.size_]),
        parent_(0),
        sibling_(0),
        children_(0),
        type_(0),
        subtype_(0) {
    for (int i = 0; i < size_; ++i)
        elements_[i] = cloneMe.elements_[i];
}

inline NTreeBag::~NTreeBag() {
    NTreeBag* tmp;
    while (children_) {
        tmp = children_;
        children_ = children_->sibling_;
        delete tmp;
    }
    delete[] elements_;
}

inline int NTreeBag::size() const {
    return size_;
}

inline int NTreeBag::element(int which) const {
    return elements_[which];
}

inline int NTreeBag::index() const {
    return index_;
}

inline int NTreeBag::type() const {
    return type_;
}

inline int NTreeBag::subtype() const {
    return subtype_;
}

inline const NTreeBag* NTreeBag::parent() const {
    return parent_;
}

inline const NTreeBag* NTreeBag::children() const {
    return children_;
}

inline const NTreeBag* NTreeBag::sibling() const {
    return sibling_;
}

inline bool NTreeBag::isLeaf() const {
    return (! children_);
}

inline void NTreeBag::insertChild(NTreeBag* child) {
    child->parent_ = this;
    child->sibling_ = children_;
    children_ = child;
}

inline void NTreeBag::swapContents(NTreeBag& other) {
    int s = size_; size_ = other.size_; other.size_ = s;
    int* e = elements_; elements_ = other.elements_; other.elements_ = e;
}

// Inline functions for NTreeDecomposition

inline NTreeDecomposition::~NTreeDecomposition() {
    delete root_;
}

inline int NTreeDecomposition::width() const {
    return width_;
}

inline int NTreeDecomposition::size() const {
    return size_;
}

inline const NTreeBag* NTreeDecomposition::root() const {
    return root_;
}

inline const NTreeBag* NTreeDecomposition::firstPrefix() const {
    return root_;
}

inline void NTreeDecomposition::reindex() {
    size_ = 0;
    for (const NTreeBag* b = first(); b; b = b->next())
        const_cast<NTreeBag*>(b)->index_ = size_++;
}

// Inline functions for NTreeDecomposition::Graph

inline NTreeDecomposition::Graph::Graph(int order) :
        order_(order), adj_(new bool*[order]) {
    int i, j;
    for (i = 0; i < order; ++i) {
        adj_[i] = new bool[order];
        for (j = 0; j < order; ++j)
            adj_[i][j] = false;
    }
}

inline NTreeDecomposition::Graph::~Graph() {
    for (int i = 0; i < order_; ++i)
        delete[] adj_[i];
    delete[] adj_;
}

} // namespace regina

#endif

