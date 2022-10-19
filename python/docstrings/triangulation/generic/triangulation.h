/*
  This file contains docstrings for use in the Python bindings.
  Do not edit! They were automatically extracted by ../gendoc.sh.
 */

#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace regina::python::doc {


// Docstring regina::python::doc::DegreeGreaterThan
static const char *DegreeGreaterThan =
R"doc(A function object used for sorting faces of triangulations by
decreasing degree. This can (for instance) be used with std::sort().

The template argument *dim* refers to the dimension of the overall
triangluation(s) with which you are working. The template argument
*subdim* refers to the dimension of the faces that you are sorting.
So, for instance, to sort edges of a 3-manifold triangulation by
decreasing edge degree, you would use DegreeGreaterThan<3, 1>.

A single instance of this class works with faces of a single fixed
triangulation (which is passed to the class constructor).

An object of this class behaves like a reference: it is lightweight
and can be copy-constructed cheaply, but it does not support
assignments or swaps.

Precondition:
    *dim* is one of Regina's standard dimensions.

Precondition:
    *subdim* is between 0 and *dim*-1 inclusive.)doc";

// Docstring regina::python::doc::DegreeLessThan
static const char *DegreeLessThan =
R"doc(A function object used for sorting faces of triangulations by
increasing degree. This can (for instance) be used with std::sort().

The template argument *dim* refers to the dimension of the overall
triangluation(s) with which you are working. The template argument
*subdim* refers to the dimension of the faces that you are sorting.
So, for instance, to sort edges of a 3-manifold triangulation by
increasing edge degree, you would use DegreeLessThan<3, 1>.

A single instance of this class works with faces of a single fixed
triangulation (which is passed to the class constructor).

An object of this class behaves like a reference: it is lightweight
and can be copy-constructed cheaply, but it does not support
assignments or swaps.

Precondition:
    *dim* is one of Regina's standard dimensions.

Precondition:
    *subdim* is between 0 and *dim*-1 inclusive.)doc";

// Docstring regina::python::doc::Triangulation
static const char *Triangulation =
R"doc(A *dim*-dimensional triangulation, built by gluing together
*dim*-dimensional simplices along their (*dim*-1)-dimensional facets.
Typically (but not necessarily) such triangulations are used to
represent *dim*-manifolds.

Such triangulations are not the same as pure simplicial complexes, for
two reasons:

* The only identifications that the user can explicitly specify are
  gluings between *dim*-dimensional simplices along their
  (*dim*-1)-dimensional facets. All other identifications between
  *k*-faces (for any *k*) are simply consequences of these
  (*dim*-1)-dimensional gluings. In contrast, a simplicial complex
  allows explicit gluings between faces of any dimension.

* There is no requirement for a *k*-face to have (*k*+1) distinct
  vertices (so, for example, edges may be loops). Many distinct
  *k*-faces of a top-dimensional simplex may be identified together as
  a consequence of the (*dim*-1)-dimensional gluings, and indeed we
  are even allowed to glue together two distinct factets of the same
  *dim*-simplex. In contrast, a simplicial complex does not allow any
  of these situations.

Amongst other things, this definition is general enough to capture any
reasonable definition of a *dim*-manifold triangulation. However,
there is no requirement that a triangulation must actually represent a
manifold (and indeed, testing this condition is undecidable for
sufficiently large *dim*).

You can construct a triangulation from scratch using routines such as
newSimplex() and Simplex<dim>::join(). There are also routines for
exporting and importing triangulations in bulk, such as isoSig() and
fromIsoSig() (which use *isomorphism signatures*), or
dumpConstruction() and fromGluings() (which use C++ code).

In additional to top-dimensional simplices, this class also tracks:

* connected components of the triangulation, as represented by the
  class Component<dim>;

* boundary components of the triangulation, as represented by the
  class BoundaryComponent<dim>;

* lower-dimensional faces of the triangulation, as represented by the
  classes Face<dim, subdim> for *subdim* = 0,...,(*dim*-1).

Such objects are temporary: whenever the triangulation changes, they
will be deleted and rebuilt, and any pointers to them will become
invalid. Likewise, if the triangulation is deleted then all component
objects will be deleted alongside it.

Since Regina 7.0, this is no longer a "packet type" that can be
inserted directly into the packet tree. Instead a Triangulation is now
a standalone mathematatical object, which makes it slimmer and faster
for ad-hoc use. The consequences of this are:

* If you create your own Triangulation, it will not have any of the
  usual packet infrastructure. You cannot add it into the packet tree,
  and it will not support a label, tags, child/parent packets, and/or
  event listeners.

* To include a Triangulation in the packet tree, you must create a new
  PacketOf<Triangulation>. This *is* a packet type, and supports
  labels, tags, child/parent packets, and event listeners. It derives
  from Triangulation, and so inherits the full Triangulation
  interface.

* If you are adding new functions to this class that edit the
  triangulation, you must still remember to create a ChangeEventSpan.
  This will ensure that, if the triangulation is being managed by a
  PacketOf<Triangulation>, then the appropriate packet change events
  will be fired. All other events (aside from packetToBeChanged() and
  packetWasChanged() are managed directly by the
  PacketOf<Triangulation> wrapper class.

This class implements C++ move semantics and adheres to the C++
Swappable requirement. It is designed to avoid deep copies wherever
possible, even when passing or returning objects by value.

For Regina's standard dimensions, this template is specialised and
offers *much* more functionality. In order to use these specialised
classes, you will need to include the corresponding headers (e.g.,
triangulation/dim2.h for *dim* = 2, or triangulation/dim3.h for *dim*
= 3).

Python:
    Python does not support templates. Instead this class can be used
    by appending the dimension as a suffix (e.g., Triangulation2 and
    Triangulation3 for dimensions 2 and 3).

Template parameter ``dim``:
    the dimension of the underlying triangulation. This must be
    between 2 and 15 inclusive.

\headerfile triangulation/generic.h)doc";

namespace DegreeGreaterThan_ {

// Docstring regina::python::doc::DegreeGreaterThan_::__call
static const char *__call =
R"doc(Compares the degrees of the *subdim*-dimensional faces at the given
indices within the working triangulation. The triangulation that is
used will be the one that was passed to the class constructor.

Precondition:
    Both *a* and *b* are non-negative, and are strictly less than the
    total number of *subdim*-dimensional faces in the triangulation.

Parameter ``a``:
    the index of the first *subdim*-dimensional face within the
    triangulation.

Parameter ``b``:
    the index of the second *subdim*-dimensional face within the
    triangulation.

Returns:
    ``True`` if and only if face *a* has greater degree than face *b*
    within the given triangulation.)doc";

// Docstring regina::python::doc::DegreeGreaterThan_::__copy
static const char *__copy = R"doc(Creates a new clone of the given function object.)doc";

// Docstring regina::python::doc::DegreeGreaterThan_::__init
static const char *__init =
R"doc(Constructions a function object for working with faces of the given
triangulation.

Parameter ``tri``:
    the triangulation with which we are working.)doc";

}

namespace DegreeLessThan_ {

// Docstring regina::python::doc::DegreeLessThan_::__call
static const char *__call =
R"doc(Compares the degrees of the *subdim*-dimensional faces at the given
indices within the working triangulation. The triangulation that is
used will be the one that was passed to the class constructor.

Precondition:
    Both *a* and *b* are non-negative, and are strictly less than the
    total number of *subdim*-dimensional faces in the triangulation.

Parameter ``a``:
    the index of the first *subdim*-dimensional face within the
    triangulation.

Parameter ``b``:
    the index of the second *subdim*-dimensional face within the
    triangulation.

Returns:
    ``True`` if and only if face *a* has smaller degree than face *b*
    within the given triangulation.)doc";

// Docstring regina::python::doc::DegreeLessThan_::__copy
static const char *__copy = R"doc(Creates a new clone of the given function object.)doc";

// Docstring regina::python::doc::DegreeLessThan_::__init
static const char *__init =
R"doc(Constructions a function object for working with faces of the given
triangulation.

Parameter ``tri``:
    the triangulation with which we are working.)doc";

}

namespace Triangulation_ {

// Docstring regina::python::doc::Triangulation_::__copy
static const char *__copy =
R"doc(Creates a new copy of the given triangulation.

This will clone any computed properties (such as homology, fundamental
group, and so on) of the given triangulation also. If you want a
"clean" copy that resets all properties to unknown, you can use the
two-argument copy constructor instead.

Parameter ``copy``:
    the triangulation to copy.)doc";

// Docstring regina::python::doc::Triangulation_::__default
static const char *__default =
R"doc(Default constructor.

Creates an empty triangulation.)doc";

// Docstring regina::python::doc::Triangulation_::__init
static const char *__init =
R"doc(Creates a new copy of the given triangulation, with the option of
whether or not to clone its computed properties also.

Parameter ``copy``:
    the triangulation to copy.

Parameter ``cloneProps``:
    ``True`` if this should also clone any computed properties of the
    given triangulation (such as homology, fundamental group, and so
    on), or ``False`` if the new triangulation should have all
    properties marked as unknown.)doc";

// Docstring regina::python::doc::Triangulation_::swap
static const char *swap =
R"doc(Swaps the contents of this and the given triangulation.

All top-dimensional simplices that belong to this triangulation will
be moved to *other*, and all top-dimensional simplices that belong to
*other* will be moved to this triangulation. Likewise, all skeletal
objects (such as lower-dimensional faces, components, and boundary
components) and all cached properties will be swapped.

In particular, any pointers or references to Simplex<dim> and/or
Face<dim, subdim> objects will remain valid.

This routine will behave correctly if *other* is in fact this
triangulation.

.. note::
    This swap function is *not* marked ``noexcept``, since it fires
    change events on both triangulations which may in turn call
    arbitrary code via any registered packet listeners.

Parameter ``other``:
    the triangulation whose contents should be swapped with this.)doc";

}

} // namespace regina::python::doc

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif
