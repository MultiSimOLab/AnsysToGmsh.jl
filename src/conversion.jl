# CDBMesh -> GmshMesh: element conversion, facet/edge extraction and
# classification of mesh entities by component membership, so that each
# ANSYS component becomes a Gmsh physical group (and therefore a Gridap tag).

"""
    MshEntity

A discrete Gmsh entity: dimension, tag, physical group tags, the nodes
classified on it (final node tags) and its elements (one element type).
"""
struct MshEntity
  dim::Int
  tag::Int
  physicals::Vector{Int}
  nodes::Vector{Int}
  etype::Int
  elements::Vector{Vector{Int}}
end

"""
    GmshMesh

Mesh ready to be written as a Gmsh 4.1 file. `coords` is `3 × nnodes` in the
final node numbering (1:nnodes). `physical_names` holds `(dim, tag, name)`.
"""
struct GmshMesh
  coords::Matrix{Float64}
  entities::Vector{MshEntity}
  physical_names::Vector{Tuple{Int,Int,String}}
end

function Base.show(io::IO, ::MIME"text/plain", m::GmshMesh)
  println(io, "GmshMesh: ", size(m.coords, 2), " nodes")
  for d in 3:-1:0
    es = filter(e -> e.dim == d, m.entities)
    isempty(es) && continue
    nel = sum(e -> length(e.elements), es)
    types = unique(e.etype for e in es if !isempty(e.elements))
    tn = join([GMSH_NAMES[t] for t in types], ", ")
    println(io, "  dim $d: ", length(es), " entities, ", nel, " elements", isempty(tn) ? "" : " ($tn)")
  end
  for d in 3:-1:0
    names = [n for (dd, _, n) in m.physical_names if dd == d]
    isempty(names) || println(io, "  physical groups dim $d: ", join(names, ", "))
  end
end

# intersection of two sorted vectors
function sorted_intersect(a::Vector{Int}, b::Vector{Int})
  out = Int[]
  i = j = 1
  while i <= length(a) && j <= length(b)
    if a[i] == b[j]
      push!(out, a[i]); i += 1; j += 1
    elseif a[i] < b[j]
      i += 1
    else
      j += 1
    end
  end
  out
end

const EMPTY_SIG = Int[]

"""
    convert_cdb(cdb::CDBMesh; kwargs...) -> GmshMesh

Build a Gmsh mesh (with physical groups) from the contents of a CDB file.

Keyword arguments
- `order = :auto`: keep quadratic Lagrangian elements (SOLID187 → tet10,
  triangular PLANE183 → tri6) and linearize serendipity elements (SOLID186
  bricks → hex8, PLANE183 quads → quad4). `order = 1` linearizes everything.
- `boundary = "boundary"`: name of a physical group with the whole boundary
  (faces, edges and vertices). `nothing` to skip it.
- `domain = "domain"`: name of a physical group with all cells. `nothing` to skip.
- `materials = "material_"`: prefix of one cell group per ANSYS material id.
  `nothing` to skip.
- `interior_faces = true`: also tag interior faces whose nodes all belong to
  a node component (interfaces between bodies). With `false` only boundary
  faces are tagged.
- `split = false`: when the mesh mixes shapes (e.g. Workbench "Hex Dominant":
  bricks, prisms, pyramids and tets), Gridap cannot load it. With
  `split = true` every cell is turned into linear tetrahedra (triangles in
  2-D) in a conforming way; non-simplex cells get an extra node at their
  centroid. Homogeneous meshes are left untouched unless `split = :always`.
- `verbose = true`: print a summary of the conversion.

Node components (CMBLOCK ...,NODE) become physical groups on vertices,
edges and faces: a face (edge) is in the group when all its corner nodes are
in the component. Element components (CMBLOCK ...,ELEM) become cell groups.
"""
function convert_cdb(cdb::CDBMesh; order=:auto, boundary="boundary", domain="domain",
                     materials="material_", interior_faces::Bool=true, split=false,
                     verbose::Bool=true)

  # ------------------------------------------------------------------
  # 1. Cells
  # ------------------------------------------------------------------
  cells = Vector{Int}[]          # ANSYS node numbers in Gmsh order
  ctypes = Int[]
  cmat = Int[]
  celid = Int[]
  skipped = Dict{Int,Int}()
  nlinear = 0
  for e in cdb.elements
    a = get(cdb.etypes, e.etype, 0)
    r = a == 0 ? nothing : ansys_to_gmsh(a, e.nodes; order)
    if r === nothing
      skipped[a] = get(skipped, a, 0) + 1
      continue
    end
    t, n, note = r
    note == :linearized && (nlinear += 1)
    push!(cells, n); push!(ctypes, t); push!(cmat, e.mat); push!(celid, e.id)
  end
  isempty(cells) && error("No supported elements found. Element types in file: " *
                          join(sort!(unique(values(cdb.etypes))), ", ") *
                          ". Supported: " * join(supported_ansys_elements(), ", "))
  D = maximum(GMSH_DIM[t] for t in ctypes)
  keep = [GMSH_DIM[t] == D for t in ctypes]
  nlower = count(!, keep)
  cells = cells[keep]; ctypes = ctypes[keep]; cmat = cmat[keep]; celid = celid[keep]

  node_col = Dict{Int,Int}()
  for (k, id) in enumerate(cdb.node_ids)
    node_col[id] = k
  end
  for c in cells, n in c
    haskey(node_col, n) || error("Element uses node $n, which is not in the NBLOCK")
  end
  X = cdb.coords

  types = unique(ctypes)
  counts = join(["$(count(==(t), ctypes)) $(GMSH_NAMES[t])" for t in types], ", ")
  needsplit = length(types) > 1 || !haskey(CELL_FACETS, types[1])
  nsplit = 0
  if needsplit && split === false
    error("The mesh mixes element shapes or uses prisms/pyramids ($counts). " *
          "Gridap/HyperFEM needs a single simplex or hypercube type per mesh. " *
          "Use `split=true` to turn every cell into linear simplices, or remesh.")
  end
  if (split === true && needsplit) || split === :always
    cells, ctypes, cmat, celid, X, nsplit = split_into_simplices(cells, ctypes, cmat, celid, X, node_col)
    verbose && @info "Split $counts into $(length(cells)) linear $(D == 3 ? "tetrahedra" : "triangles")"
  end
  T = ctypes[1]
  FT = FACET_TYPE[T]
  nverts = Dict(GMSH_TET4 => 4, GMSH_TET10 => 4, GMSH_HEX8 => 8,
                GMSH_TRI3 => 3, GMSH_TRI6 => 3, GMSH_QUAD4 => 4)[T]

  corner = Set{Int}()
  for c in cells
    for k in 1:nverts
      push!(corner, c[k])
    end
  end

  midmap = Dict{Tuple{Int,Int},Int}()
  if isquadratic(T)
    for c in cells, (p, m) in CELL_MIDNODES[T]
      a, b = c[p[1]], c[p[2]]
      midmap[minmax(a, b)] = c[m]
    end
  end
  mid(a, b) = midmap[minmax(a, b)]

  # ------------------------------------------------------------------
  # 2. Node component signatures (component indices are 1-based;
  #    index 0 is reserved for the automatic boundary group)
  # ------------------------------------------------------------------
  ncomp_names = first.(cdb.node_components)
  use_boundary = boundary !== nothing && !(boundary in ncomp_names)
  if boundary !== nothing && !use_boundary
    @warn "A node component is already called '$boundary'; the automatic boundary group is not added."
  end
  node_sig = Dict{Int,Vector{Int}}()
  for (k, (name, ids)) in enumerate(cdb.node_components)
    for id in ids
      haskey(node_col, id) || continue
      s = get!(node_sig, id, Int[])
      (isempty(s) || s[end] != k) && push!(s, k)
    end
  end
  nsig(n) = get(node_sig, n, EMPTY_SIG)

  # ------------------------------------------------------------------
  # 3. Facets (faces in 3-D, edges in 2-D)
  # ------------------------------------------------------------------
  facet_id = Dict{Vector{Int},Int}()
  facet_nodes = Vector{Int}[]
  facet_count = Int[]
  for c in cells
    for lf in CELL_FACETS[T]
      cn = [c[i] for i in lf]
      key = sort(cn)
      f = get(facet_id, key, 0)
      if f == 0
        if FT == GMSH_TRI6
          a, b, d = cn
          push!(facet_nodes, [a, b, d, mid(a, b), mid(b, d), mid(d, a)])
        elseif FT == GMSH_LINE3
          push!(facet_nodes, [cn[1], cn[2], mid(cn[1], cn[2])])
        else
          push!(facet_nodes, cn)
        end
        push!(facet_count, 1)
        facet_id[key] = length(facet_nodes)
      else
        facet_count[f] += 1
      end
    end
  end
  nfc = length(CELL_FACETS[T][1])            # corners per facet
  isbnd(f) = facet_count[f] == 1

  bnodes = Set{Int}()
  for f in eachindex(facet_nodes)
    isbnd(f) && union!(bnodes, facet_nodes[f])
  end

  function corner_sig(nodes)
    s = nsig(nodes[1])
    for k in 2:length(nodes)
      isempty(s) && break
      s = sorted_intersect(s, nsig(nodes[k]))
    end
    s
  end
  withbnd(s, onb) = (use_boundary && onb) ? vcat(0, s) : s

  facet_sig = Vector{Vector{Int}}(undef, length(facet_nodes))
  for f in eachindex(facet_nodes)
    b = isbnd(f)
    s = (b || interior_faces) ? corner_sig(view(facet_nodes[f], 1:nfc)) : Int[]
    facet_sig[f] = withbnd(s, b)
  end

  # ------------------------------------------------------------------
  # 4. Edges (3-D only): edges of boundary faces and of tagged faces
  # ------------------------------------------------------------------
  edge_nodes = Vector{Int}[]
  edge_sig = Vector{Vector{Int}}()
  if D == 3
    edge_id = Dict{Tuple{Int,Int},Int}()
    edge_onb = Bool[]
    for f in eachindex(facet_nodes)
      (isbnd(f) || !isempty(facet_sig[f])) || continue
      fn = facet_nodes[f]
      for k in 1:nfc
        a, b = fn[k], fn[mod1(k + 1, nfc)]
        key = minmax(a, b)
        e = get(edge_id, key, 0)
        if e == 0
          push!(edge_nodes, FT == GMSH_TRI6 ? [a, b, mid(a, b)] : [a, b])
          push!(edge_onb, isbnd(f))
          edge_id[key] = length(edge_nodes)
        else
          edge_onb[e] |= isbnd(f)
        end
      end
    end
    for e in eachindex(edge_nodes)
      push!(edge_sig, withbnd(corner_sig(view(edge_nodes[e], 1:2)), edge_onb[e]))
    end
  end

  # ------------------------------------------------------------------
  # 5. Entities: one per distinct signature and dimension
  # ------------------------------------------------------------------
  function group_by_sig(items, sigs)
    tags = Dict{Vector{Int},Int}()
    members = Vector{Vector{Int}}()
    keys_ = Vector{Vector{Int}}()
    for (i, s) in zip(items, sigs)
      isempty(s) && continue
      t = get!(tags, s) do
        push!(members, Int[]); push!(keys_, s); length(members)
      end
      push!(members[t], i)
    end
    keys_, members
  end

  cn_list = sort!(collect(corner))
  pt_sigs, pt_members = group_by_sig(cn_list, [withbnd(nsig(n), n in bnodes) for n in cn_list])
  fc_sigs, fc_members = group_by_sig(eachindex(facet_nodes), facet_sig)
  ed_sigs, ed_members = D == 3 ? group_by_sig(eachindex(edge_nodes), edge_sig) : (Vector{Int}[], Vector{Int}[])

  # cell entities: (material, element components)
  ecomp_names = first.(cdb.element_components)
  cell_ecomps = Dict{Int,Vector{Int}}()
  for (k, (_, ids)) in enumerate(cdb.element_components), id in ids
    s = get!(cell_ecomps, id, Int[])
    (isempty(s) || s[end] != k) && push!(s, k)
  end
  vol_key = Dict{Tuple{Int,Vector{Int}},Int}()
  vol_keys = Tuple{Int,Vector{Int}}[]
  vol_members = Vector{Int}[]
  for i in eachindex(cells)
    key = (cmat[i], get(cell_ecomps, celid[i], EMPTY_SIG))
    t = get!(vol_key, key) do
      push!(vol_keys, key); push!(vol_members, Int[]); length(vol_members)
    end
    push!(vol_members[t], i)
  end

  # ------------------------------------------------------------------
  # 6. Physical groups
  # ------------------------------------------------------------------
  physical_names = Tuple{Int,Int,String}[]
  ent_phys = Dict{Tuple{Int,Int},Vector{Int}}()     # (dim, entity) -> physical tags
  addphys!(d, ents, ptag) = for e in ents; push!(get!(ent_phys, (d, e), Int[]), ptag); end

  lower = D == 3 ? [(0, pt_sigs), (1, ed_sigs), (2, fc_sigs)] : [(0, pt_sigs), (1, fc_sigs)]
  comp_list = Tuple{Int,String}[]
  use_boundary && push!(comp_list, (0, boundary))
  for (k, name) in enumerate(ncomp_names)
    push!(comp_list, (k, name))
  end
  for (d, sigs) in lower
    ptag = 0
    for (k, name) in comp_list
      ents = [t for (t, s) in enumerate(sigs) if insorted(k, s)]
      isempty(ents) && continue
      ptag += 1
      push!(physical_names, (d, ptag, name))
      addphys!(d, ents, ptag)
    end
  end

  ptag = 0
  allvol = collect(eachindex(vol_keys))
  if domain !== nothing && !(domain in ecomp_names)
    ptag += 1; push!(physical_names, (D, ptag, domain)); addphys!(D, allvol, ptag)
  end
  if materials !== nothing
    for m in sort!(unique(cmat))
      ptag += 1
      push!(physical_names, (D, ptag, string(materials, m)))
      addphys!(D, [t for (t, k) in enumerate(vol_keys) if k[1] == m], ptag)
    end
  end
  for (k, name) in enumerate(ecomp_names)
    ents = [t for (t, key) in enumerate(vol_keys) if insorted(k, key[2])]
    if isempty(ents)
      verbose && @warn "Element component '$name' contains no converted cells; skipped."
      continue
    end
    ptag += 1; push!(physical_names, (D, ptag, name)); addphys!(D, ents, ptag)
  end

  # ------------------------------------------------------------------
  # 7. Node numbering: point entities first, then cell entities, so that
  #    node tags are consecutive in the order Gmsh returns them.
  # ------------------------------------------------------------------
  newtag = Dict{Int,Int}()
  ent_nodes = Dict{Tuple{Int,Int},Vector{Int}}()
  order_ids = Int[]
  for (t, ns) in enumerate(pt_members)
    lst = Int[]
    for n in ns
      haskey(newtag, n) && continue
      push!(order_ids, n); newtag[n] = length(order_ids); push!(lst, newtag[n])
    end
    ent_nodes[(0, t)] = lst
  end
  for (t, cs) in enumerate(vol_members)
    lst = Int[]
    for i in cs, n in cells[i]
      haskey(newtag, n) && continue
      push!(order_ids, n); newtag[n] = length(order_ids); push!(lst, newtag[n])
    end
    ent_nodes[(D, t)] = lst
  end
  coords = Matrix{Float64}(undef, 3, length(order_ids))
  for (k, n) in enumerate(order_ids)
    coords[:, k] = view(X, :, node_col[n])
  end
  rn(v) = [newtag[n] for n in v]

  # ------------------------------------------------------------------
  # 8. Entities
  # ------------------------------------------------------------------
  entities = MshEntity[]
  mk(d, t, etype, elems) = MshEntity(d, t, get(ent_phys, (d, t), Int[]), get(ent_nodes, (d, t), Int[]), etype, elems)
  for (t, ns) in enumerate(pt_members)
    push!(entities, mk(0, t, GMSH_POINT, [[newtag[n]] for n in ns]))
  end
  if D == 3
    ET = FACET_EDGE_TYPE[FT]
    for (t, es) in enumerate(ed_members)
      push!(entities, mk(1, t, ET, [rn(edge_nodes[e]) for e in es]))
    end
  end
  for (t, fs) in enumerate(fc_members)
    push!(entities, mk(D - 1, t, FT, [rn(facet_nodes[f]) for f in fs]))
  end
  for (t, cs) in enumerate(vol_members)
    push!(entities, mk(D, t, T, [rn(cells[i]) for i in cs]))
  end

  mesh = GmshMesh(coords, entities, physical_names)

  if verbose
    msg = IOBuffer()
    println(msg, "Converted $(length(cells)) cells ($(GMSH_NAMES[T])), $(length(order_ids)) nodes")
    nlinear > 0 && println(msg, "  $nlinear quadratic serendipity elements were made linear (midside nodes dropped)")
    nlower > 0 && println(msg, "  $nlower lower-dimensional elements were ignored")
    for (a, n) in sort!(collect(skipped))
      println(msg, "  $n elements of unsupported type $(a == 0 ? "(undefined ET)" : a) were ignored")
    end
    for (k, name) in enumerate(ncomp_names)
      any(s -> insorted(k, s), fc_sigs) && continue
      any(s -> insorted(k, s), pt_sigs) || (println(msg, "  node component '$name' has no mesh vertices"); continue)
      println(msg, "  node component '$name' only tags vertices/edges (no full face)")
    end
    show(msg, MIME("text/plain"), mesh)
    @info String(take!(msg))
  end
  mesh
end

# ---------------------------------------------------------------------------
# Conforming split of any linear cell into simplices
# ---------------------------------------------------------------------------

function tet_volume(X, node_col, t)
  a, b, c, d = (view(X, :, node_col[n]) for n in t)
  orient3(a, b, c, d)
end

# 6 × signed volume of the tet (a,b,c,d)
function orient3(a, b, c, d)
  u1, u2, u3 = b[1] - a[1], b[2] - a[2], b[3] - a[3]
  v1, v2, v3 = c[1] - a[1], c[2] - a[2], c[3] - a[3]
  w1, w2, w3 = d[1] - a[1], d[2] - a[2], d[3] - a[3]
  u1 * (v2 * w3 - v3 * w2) - u2 * (v1 * w3 - v3 * w1) + u3 * (v1 * w2 - v2 * w1)
end

tri_area(X, node_col, t) = begin
  a, b, c = (view(X, :, node_col[n]) for n in t)
  (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
end

# Triangulated boundary of a linear cell, outward oriented (node ids)
function cell_triangles(t, c)
  tris = NTuple{3,Int}[]
  for lf in SPLIT_FACES[t]
    f = [c[i] for i in lf]
    if length(f) == 3
      push!(tris, Tuple(f))
    else
      append!(tris, split_quad(f))
    end
  end
  tris
end

# smallest 6×volume of the tets (tri, x); positive = valid split
min_cone_volume(tris, P, x) = minimum(orient3(P(t[1]), P(t[3]), P(t[2]), x) for t in tris)

# Solve the 4×4 system A z = b by Gaussian elimination; nothing if singular
function solve4(A, b)
  M = hcat(A, b)
  for k in 1:4
    p = k - 1 + argmax(abs.(M[k:4, k]))
    abs(M[p, k]) < 1e-14 && return nothing
    M[k, :], M[p, :] = M[p, :], M[k, :]
    for i in k+1:4
      M[i, :] .-= (M[i, k] / M[k, k]) .* M[k, :]
    end
  end
  z = zeros(4)
  for k in 4:-1:1
    z[k] = (M[k, 5] - sum(M[k, j] * z[j] for j in k+1:4; init=0.0)) / M[k, k]
  end
  z
end

# Centre of the largest ball inside all face planes (the kernel of a
# star-shaped cell): maximise r s.t. n_i⋅x + r ≤ n_i⋅a_i. Small LP solved by
# enumerating its vertices. Returns (x, r).
function kernel_center(tris, P)
  N = Vector{Float64}[]; d = Float64[]
  for t in tris
    a, b, c = P(t[1]), P(t[2]), P(t[3])
    n = [(b[2]-a[2])*(c[3]-a[3]) - (b[3]-a[3])*(c[2]-a[2]),
         (b[3]-a[3])*(c[1]-a[1]) - (b[1]-a[1])*(c[3]-a[3]),
         (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])]
    l = sqrt(sum(abs2, n)); l == 0 && continue
    n ./= l
    push!(N, n); push!(d, n[1]*a[1] + n[2]*a[2] + n[3]*a[3])
  end
  m = length(N)
  best = (zeros(3), -Inf)
  for i in 1:m, j in i+1:m, k in j+1:m, l in k+1:m
    idx = (i, j, k, l)
    A = [N[q][r] for q in idx, r in 1:3]
    z = solve4(hcat(A, ones(4)), [d[q] for q in idx])
    z === nothing && continue
    x = z[1:3]; r = z[4]
    r <= best[2] && continue
    all(N[q][1]*x[1] + N[q][2]*x[2] + N[q][3]*x[3] + r <= d[q] + 1e-9 * (1 + abs(d[q])) for q in 1:m) || continue
    best = (x, r)
  end
  best
end

"""
Split every cell into linear simplices, conformingly. 3-D: tets are kept
(made linear); pyramids are cut into two tets along the base diagonal
through the smallest node id; other cells get one extra node inside them
(their vertex average or, if that does not give valid tets, the centre of
their kernel) and one tet per triangulated face. Quadrilateral faces are
always cut along the diagonal through their smallest node id, so
neighbouring cells agree. 2-D: quads are cut the same way.
Returns (cells, ctypes, cmat, celid, X, nnew), X extended by the new nodes
(registered in `node_col`).
"""
function split_into_simplices(cells, ctypes, cmat, celid, X, node_col)
  ncells = Vector{Int}[]; nt = Int[]; nm = Int[]; nid = Int[]
  nextid = maximum(keys(node_col)) + 1
  ncol = size(X, 2)
  newcoords = Vector{Float64}[]
  P(n) = view(X, :, node_col[n])
  for (c0, t0, m, id) in zip(cells, ctypes, cmat, celid)
    t, c = linear_corners(t0, c0)
    out = Vector{Int}[]
    if t == GMSH_TET4 || t == GMSH_TRI3
      push!(out, collect(c))
    elseif t == GMSH_QUAD4
      for tri in split_quad(c)
        push!(out, collect(tri))
      end
    elseif haskey(SPLIT_FACES, t)
      done = false
      if t == GMSH_PYR5
        # base (1,2,3,4) seen from the apex is clockwise in SPLIT_FACES order
        base = split_quad([c[1], c[4], c[3], c[2]])
        cand = [[tri[1], tri[3], tri[2], c[5]] for tri in base]
        if all(s -> orient3(P.(s)...) > 0, cand)
          append!(out, cand); done = true
        end
      end
      if !done
        tris = cell_triangles(t, c)
        x = sum(P(n) for n in c) ./ length(c)
        if min_cone_volume(tris, P, x) <= 0
          x, r = kernel_center(tris, P)
          r > 0 || error("Cannot split $(GMSH_NAMES[t]) element $id into valid tetrahedra (not star-shaped)")
        end
        push!(newcoords, x)
        cid = nextid; nextid += 1
        ncol += 1; node_col[cid] = ncol
        for tri in tris
          push!(out, [tri[1], tri[3], tri[2], cid])
        end
      end
    else
      error("Cannot split $(GMSH_NAMES[t]) cells")
    end
    for s in out
      push!(ncells, s); push!(nt, length(s) == 4 ? GMSH_TET4 : GMSH_TRI3); push!(nm, m); push!(nid, id)
    end
  end
  Xn = isempty(newcoords) ? X : hcat(X, reduce(hcat, newcoords))
  nneg = 0
  for s in ncells
    v = length(s) == 4 ? tet_volume(Xn, node_col, s) : tri_area(Xn, node_col, s)
    if v < 0
      nneg += 1
      s[1], s[2] = s[2], s[1]
    end
  end
  nneg > 0 && @warn "$nneg simplices had negative orientation and were reoriented (check the source mesh)."
  ncells, nt, nm, nid, Xn, length(newcoords)
end
