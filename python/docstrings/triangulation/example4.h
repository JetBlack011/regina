/*
  This file contains docstrings for use in the Python bindings.
  Do not edit! They were automatically extracted by ../gendoc.sh.
 */

#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace regina::python::doc {


// Docstring regina::python::doc::Example
static const char *Example =
R"doc(Offers routines for constructing a variety of sample 4-dimensional
triangulations.

This is a specialisation of the generic Example class template; see
the Example template documentation for a general overview of how the
example triangulation classes work.

This 4-dimensional specialisation offers significant extra
functionality, by providing several more hard-coded and parameterised
constructions.)doc";

namespace Example_ {

// Docstring regina::python::doc::Example_::bundleWithMonodromy
static const char *bundleWithMonodromy =
R"doc(Returns a bundle formed from a given 3-manifold and a given monodromy.

Specifically, let *M* be the given 3-manifold triangulation. This
routine builds the bundle ``M x I``, and then identifies the two
copies of *M* on the boundary according to the given homeomorphism
from *M* to itself. The homeomorphism must be expressed as a
combinatorial automorphism, which means that for a non-trivial
monodromy you may need to do some work to find a sufficiently
symmetric 3-manifold triangulation to begin with.

The resulting manifold will contain 82 pentachora for each original
tetrahedron of *M*, and will contain many internal vertices. It is
highly recommended that you call
Triangulation<4>::intelligentSimplify() afterwards if you do not need
to preserve the combinatorial structure.

Precondition:
    The given monodromy must be an isomorphism from *M* to itself;
    that is, a combinatorial automorphism.

.. warning::
    If the given 3-manifold triangulation has ideal boundary, then you
    will obtain an invalid 4-manifold triangulation as a result.

Parameter ``base``:
    the 3-manifold triangulation *M*, as described above.

Parameter ``monodromy``:
    the homeomorphism from *M* to itself, as described above.

Returns:
    the requested bundle.)doc";

// Docstring regina::python::doc::Example_::cappellShaneson
static const char *cappellShaneson =
R"doc(Returns a two-pentachoron triangulation of a Cappell-Shaneson 2-knot
complement in the 4-sphere. This triangulation is described and
analysed in "Triangulating a Cappell-Shaneson knot complement",
Budney, Burton and Hillman, Mathematical Research Letters 19 (2012),
no. 5, 1117-1126.

Returns:
    a Cappell-Shaneson 2-knot complement.)doc";

// Docstring regina::python::doc::Example_::cp2
static const char *cp2 =
R"doc(Returns a four-pentachoron triangulation of the standard complex
projective plane. This triangulation is minimal.

Under the orientation convention that we use for intersection forms,
this triangulation gives the "plain" ``CP^2`` with intersection form
[1], not the reflected ``CP^2`` with intersection form [-1].

Returns:
    the standard complex projective plane.)doc";

// Docstring regina::python::doc::Example_::fourSphere
static const char *fourSphere =
R"doc(Returns a two-pentachoron triangulation of the 4-sphere. This is
identical to calling the generic routine sphere().

Returns:
    a two-pentachoron 4-sphere.)doc";

// Docstring regina::python::doc::Example_::iBundle
static const char *iBundle =
R"doc(Returns a triangulation of the product ``M x I``, where *M* is the
given 3-manifold triangulation.

The boundary of this product will consist of two copies of *M*, both
combinatorially isomorphic to the original triangulation. If *n* is
the number of tetrahedra in *M*, then the first copy of *M* on the
boundary is obtained by mapping vertices 0,1,2,3 of tetrahedron *i* of
*M* to vertices 0,1,2,3 of pentachoron *i*, and the second copy is
obtained by mapping vertices 0,1,2,3 of tetrahedron *i* of *M* to
vertices 0,1,2,3 of pentachoron *n*+i.

The product itself will contain 82 pentachora for each original
tetrahedron of *M*, and will contain many internal vertices. It is
highly recommended that you call
Triangulation<4>::intelligentSimplify() afterwards if you do not need
to preserve the combinatorial structure.

.. warning::
    If the given 3-manifold triangulation has ideal boundary, then you
    will obtain an invalid 4-manifold triangulation as a result.

Parameter ``base``:
    the 3-manifold triangulation *M*, as described above.

Returns:
    the product ``M x I``.)doc";

// Docstring regina::python::doc::Example_::k3
static const char *k3 =
R"doc(Returns a triangulation of the standard K3 surface.

Be warned: this triangulation is extremely large.

Returns:
    the K3 surface.)doc";

// Docstring regina::python::doc::Example_::rp4
static const char *rp4 =
R"doc(Returns a four-pentachoron triangulation of real projective 4-space.

Returns:
    real projective 4-space.)doc";

// Docstring regina::python::doc::Example_::s1Bundle
static const char *s1Bundle =
R"doc(Returns a triangulation of the product ``M x S1``, where *M* is the
given 3-manifold triangulation. This simply calls iBundle() and then
glues together the two copies of *M* on the boundary.

The product will contain 82 pentachora for each original tetrahedron
of *M*, and will contain many internal vertices. It is highly
recommended that you call Triangulation<4>::intelligentSimplify()
afterwards if you do not need to preserve the combinatorial structure.

.. warning::
    If the given 3-manifold triangulation has ideal boundary, then you
    will obtain an invalid 4-manifold triangulation as a result.

Parameter ``base``:
    the 3-manifold triangulation *M*, as described above.

Returns:
    the product ``M x S1``.)doc";

// Docstring regina::python::doc::Example_::s2xs2
static const char *s2xs2 =
R"doc(Returns a six-pentachoron triangulation of the standard product ``S^2
x S^2``. This triangulation is minimal.

Returns:
    the standard product of two 2-spheres.)doc";

// Docstring regina::python::doc::Example_::s2xs2Twisted
static const char *s2xs2Twisted =
R"doc(Returns a six-pentachoron triangulation of the twisted product ``S^2
x~ S^2``. This manifold is diffeomorphic to ``CP^2 # -CP^2``, where
``-CP^2`` denotes ``CP^2`` with its orientation reversed. This
triangulation is minimal.

Returns:
    the twisted product of two 2-spheres.)doc";

// Docstring regina::python::doc::Example_::s3xs1
static const char *s3xs1 =
R"doc(Returns a two-pentachoron triangulation of the product space ``S^3 x
S^1``. This is identical to calling the generic routine
sphereBundle().

Returns:
    the product ``S^3 x S^1``.)doc";

// Docstring regina::python::doc::Example_::s3xs1Twisted
static const char *s3xs1Twisted =
R"doc(Returns a two-pentachoron triangulation of the twisted product space
``S^3 x~ S^1``. This is identical to calling the generic routine
twistedSphereBundle().

Returns:
    the twisted product ``S^3 x~ S^1``.)doc";

// Docstring regina::python::doc::Example_::simplicialFourSphere
static const char *simplicialFourSphere =
R"doc(Returns the standard six-pentachoron triangulation of the 4-sphere as
the boundary of a 5-simplex. This is identical to calling the generic
routine simplicialSphere().

Returns:
    the standard simplicial 4-sphere.)doc";

}

} // namespace regina::python::doc

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

