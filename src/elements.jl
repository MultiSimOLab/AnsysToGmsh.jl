# ANSYS element -> Gmsh element conversion.
#
# Gmsh element types used here
const GMSH_LINE2 = 1
const GMSH_TRI3 = 2
const GMSH_QUAD4 = 3
const GMSH_TET4 = 4
const GMSH_HEX8 = 5
const GMSH_PRISM6 = 6
const GMSH_PYR5 = 7
const GMSH_LINE3 = 8
const GMSH_TRI6 = 9
const GMSH_TET10 = 11
const GMSH_POINT = 15

const GMSH_NAMES = Dict(
  GMSH_LINE2 => "line2", GMSH_TRI3 => "tri3", GMSH_QUAD4 => "quad4",
  GMSH_TET4 => "tet4", GMSH_HEX8 => "hex8", GMSH_PRISM6 => "prism6",
  GMSH_PYR5 => "pyramid5", GMSH_LINE3 => "line3", GMSH_TRI6 => "tri6",
  GMSH_TET10 => "tet10", GMSH_POINT => "point")

const GMSH_DIM = Dict(
  GMSH_POINT => 0, GMSH_LINE2 => 1, GMSH_LINE3 => 1,
  GMSH_TRI3 => 2, GMSH_QUAD4 => 2, GMSH_TRI6 => 2,
  GMSH_TET4 => 3, GMSH_HEX8 => 3, GMSH_PRISM6 => 3, GMSH_PYR5 => 3, GMSH_TET10 => 3)

# Families of ANSYS elements that share the same node layout.
# 8-node bricks: I J K L (bottom) M N O P (top)
const ANSYS_HEX8 = (5, 45, 65, 69, 70, 96, 164, 185, 278)
# 20-node bricks: I..P corners, Q R S T, U V W X, Y Z A B midside nodes
const ANSYS_HEX20 = (90, 95, 117, 122, 186, 226, 231, 236, 279)
# 10-node tets: I J K L, M(IJ) N(JK) O(KI) P(IL) Q(JL) R(KL)
const ANSYS_TET10 = (87, 92, 98, 123, 187, 227, 232, 237)
# 4-node tets
const ANSYS_TET4 = (72, 285)
# 4-node quads (2-D): I J K L
const ANSYS_QUAD4 = (13, 42, 55, 182)
# 8-node quads (2-D): I J K L, M(IJ) N(JK) O(KL) P(LI)
const ANSYS_QUAD8 = (53, 77, 82, 121, 183, 223, 230, 233)
# 6-node triangles (2-D): I J K, L(IJ) M(JK) N(KI)
const ANSYS_TRI6 = (2, 35)

"""
    supported_ansys_elements()

ANSYS element numbers that `AnsysToGmsh` knows how to convert.
"""
supported_ansys_elements() =
  sort!(collect((ANSYS_HEX8..., ANSYS_HEX20..., ANSYS_TET10..., ANSYS_TET4...,
                 ANSYS_QUAD4..., ANSYS_QUAD8..., ANSYS_TRI6...)))

allsame(v) = all(==(v[1]), v)

"""
    ansys_to_gmsh(ansys_type, nodes; order=:auto) -> (gmsh_type, gmsh_nodes, note)

Convert one ANSYS element (ANSYS node order) to a Gmsh element. Returns
`nothing` when the element type is not supported. `note` is a symbol that
records lossy conversions (`:linearized`) or `:none`.

`order = :auto` keeps quadratic Lagrangian elements (tet10, tri6) and turns
serendipity elements (hex20, quad8) into linear ones, since Gridap only
handles Lagrangian elements. `order = 1` makes every element linear.
"""
function ansys_to_gmsh(a::Integer, n::AbstractVector{<:Integer}; order=:auto)
  order in (:auto, 1) || throw(ArgumentError("order must be :auto or 1, got $order"))
  lin = order == 1
  if a in ANSYS_HEX8 || a in ANSYS_HEX20
    quad = a in ANSYS_HEX20
    length(n) >= (quad ? 20 : 8) || error("ANSYS element $a needs $(quad ? 20 : 8) nodes, got $(length(n))")
    c = n[1:8]
    istet = c[3] == c[4] && allsame(c[5:8])
    if istet
      if quad && !lin
        # degenerate 20-node brick = 10-node tet
        return (GMSH_TET10, n[[1, 2, 3, 5, 9, 10, 12, 17, 19, 18]], :none)
      end
      return (GMSH_TET4, c[[1, 2, 3, 5]], quad ? :linearized : :none)
    elseif allsame(c[5:8])
      return (GMSH_PYR5, c[1:5], quad ? :linearized : :none)
    elseif c[3] == c[4] && c[7] == c[8]
      return (GMSH_PRISM6, c[[1, 2, 3, 5, 6, 7]], quad ? :linearized : :none)
    elseif length(unique(c)) == 8
      return (GMSH_HEX8, c, quad ? :linearized : :none)
    else
      error("Unsupported degenerate form of ANSYS element $a with corner nodes $c")
    end
  elseif a in ANSYS_TET10
    length(n) >= 10 || error("ANSYS element $a needs 10 nodes, got $(length(n))")
    lin && return (GMSH_TET4, n[1:4], :linearized)
    return (GMSH_TET10, n[[1, 2, 3, 4, 5, 6, 7, 8, 10, 9]], :none)
  elseif a in ANSYS_TET4
    return (GMSH_TET4, n[[1, 2, 3, 4]], :none)
  elseif a in ANSYS_QUAD4 || a in ANSYS_QUAD8
    quad = a in ANSYS_QUAD8
    c = n[1:4]
    if c[3] == c[4]
      if quad && !lin
        return (GMSH_TRI6, n[[1, 2, 3, 5, 6, 8]], :none)
      end
      return (GMSH_TRI3, c[1:3], quad ? :linearized : :none)
    end
    length(unique(c)) == 4 || error("Unsupported degenerate form of ANSYS element $a with corner nodes $c")
    return (GMSH_QUAD4, c, quad ? :linearized : :none)
  elseif a in ANSYS_TRI6
    lin && return (GMSH_TRI3, n[1:3], :linearized)
    return (GMSH_TRI6, n[1:6], :none)
  end
  nothing
end

# ---------------------------------------------------------------------------
# Local topology of the Gmsh elements (1-based local node numbers)
# ---------------------------------------------------------------------------

# corner nodes of the faces (3-D) / edges (2-D) of each cell type
const CELL_FACETS = Dict(
  GMSH_TET4 => [(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)],
  GMSH_TET10 => [(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)],
  GMSH_HEX8 => [(1, 4, 3, 2), (1, 2, 6, 5), (1, 5, 8, 4), (2, 3, 7, 6), (3, 4, 8, 7), (5, 6, 7, 8)],
  GMSH_TRI3 => [(1, 2), (2, 3), (3, 1)],
  GMSH_TRI6 => [(1, 2), (2, 3), (3, 1)],
  GMSH_QUAD4 => [(1, 2), (2, 3), (3, 4), (4, 1)])

# local edge (corner pair) -> local midside node, quadratic cells only
const CELL_MIDNODES = Dict(
  GMSH_TET10 => [(1, 2) => 5, (2, 3) => 6, (3, 1) => 7, (4, 1) => 8, (4, 3) => 9, (4, 2) => 10],
  GMSH_TRI6 => [(1, 2) => 4, (2, 3) => 5, (3, 1) => 6])

# facet element type for each cell type
const FACET_TYPE = Dict(
  GMSH_TET4 => GMSH_TRI3, GMSH_TET10 => GMSH_TRI6, GMSH_HEX8 => GMSH_QUAD4,
  GMSH_TRI3 => GMSH_LINE2, GMSH_TRI6 => GMSH_LINE3, GMSH_QUAD4 => GMSH_LINE2)

# edge element type of a 3-D facet type
const FACET_EDGE_TYPE = Dict(GMSH_TRI3 => GMSH_LINE2, GMSH_TRI6 => GMSH_LINE3, GMSH_QUAD4 => GMSH_LINE2)

isquadratic(t) = t in (GMSH_TET10, GMSH_TRI6, GMSH_LINE3)

# ---------------------------------------------------------------------------
# Splitting into simplices (for meshes that mix shapes, e.g. Hex Dominant)
# ---------------------------------------------------------------------------

# faces of linear 3-D cells, oriented outwards (right-hand rule)
const SPLIT_FACES = Dict(
  GMSH_HEX8 => [(1, 4, 3, 2), (1, 2, 6, 5), (1, 5, 8, 4), (2, 3, 7, 6), (3, 4, 8, 7), (5, 6, 7, 8)],
  GMSH_PRISM6 => [(1, 3, 2), (4, 5, 6), (1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6)],
  GMSH_PYR5 => [(1, 4, 3, 2), (1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 1, 5)])

# Linear corner nodes of a converted cell
function linear_corners(t, n)
  t == GMSH_TET10 && return (GMSH_TET4, n[1:4])
  t == GMSH_TRI6 && return (GMSH_TRI3, n[1:3])
  (t, n)
end

# Split a quadrilateral (cyclic order) along the diagonal through its
# smallest node id. Neighbouring cells pick the same diagonal, so the
# resulting triangulation is conforming.
function split_quad(q)
  k = argmin(q)
  a, b, c, d = q[k], q[mod1(k + 1, 4)], q[mod1(k + 2, 4)], q[mod1(k + 3, 4)]
  ((a, b, c), (a, c, d))
end
