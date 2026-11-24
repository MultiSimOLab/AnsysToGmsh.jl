# Reader for ANSYS Mechanical APDL archive files (.cdb), as written by CDWRITE
# or exported from Workbench. Only the blocks needed to build a mesh are read:
#   ET       element type id -> ANSYS element number (e.g. 1 -> 187)
#   NBLOCK   node coordinates
#   EBLOCK   element connectivity (SOLID and non-SOLID formats)
#   CMBLOCK  node and element components (named selections)
# Everything else (materials, loads, keyopts, ...) is ignored.

"""
    CDBElement

One element read from an EBLOCK: ANSYS element number `id`, element type id
`etype` (local id, see `CDBMesh.etypes`), material `mat`, real constant `real`,
section `secnum` and the ANSYS node numbers in ANSYS order.
"""
struct CDBElement
  id::Int
  etype::Int
  mat::Int
  real::Int
  secnum::Int
  nodes::Vector{Int}
end

"""
    CDBMesh

Raw contents of a CDB file.

- `node_ids`, `coords`: ANSYS node numbers and their coordinates (`3 × n`).
- `elements`: all elements of the EBLOCKs (supported or not).
- `etypes`: element type id → ANSYS element number (`1 => 187`).
- `node_components`, `element_components`: component name → ANSYS ids,
  in the order they appear in the file.
"""
struct CDBMesh
  node_ids::Vector{Int}
  coords::Matrix{Float64}
  elements::Vector{CDBElement}
  etypes::Dict{Int,Int}
  node_components::Vector{Pair{String,Vector{Int}}}
  element_components::Vector{Pair{String,Vector{Int}}}
end

function Base.show(io::IO, ::MIME"text/plain", m::CDBMesh)
  println(io, "CDBMesh: ", length(m.node_ids), " nodes, ", length(m.elements), " elements")
  counts = Dict{Int,Int}()
  for e in m.elements
    a = get(m.etypes, e.etype, 0)
    counts[a] = get(counts, a, 0) + 1
  end
  for a in sort!(collect(keys(counts)))
    name = a == 0 ? "undefined type" : "ANSYS element $a"
    println(io, "  ", rpad(name, 22), counts[a])
  end
  isempty(m.node_components) || println(io, "  node components:    ", join(first.(m.node_components), ", "))
  isempty(m.element_components) || print(io, "  element components: ", join(first.(m.element_components), ", "))
end

# ---------------------------------------------------------------------------
# Fixed-width (Fortran format) helpers
# ---------------------------------------------------------------------------

# "(3i9,6e21.13e3)" -> [(:i,9),(:i,9),(:i,9),(:e,21),...]
function parse_fortran_format(line::AbstractString)
  s = strip(line)
  s = strip(s, ['(', ')'])
  fields = Tuple{Symbol,Int}[]
  for item in split(s, ',')
    item = strip(item)
    isempty(item) && continue
    m = match(r"^(\d*)\s*([iIeEfFgGdDaA])\s*(\d+)", item)
    m === nothing && error("Unsupported Fortran format item '$item' in '$line'")
    rep = isempty(m.captures[1]) ? 1 : parse(Int, m.captures[1])
    kind = lowercase(m.captures[2]) == "i" ? :i : :r
    w = parse(Int, m.captures[3])
    for _ in 1:rep
      push!(fields, (kind, w))
    end
  end
  isempty(fields) && error("Empty Fortran format '$line'")
  fields
end

# Substrings of `line` according to the field widths. Stops at the end of the line.
function fixed_width_split(line::AbstractString, fmt)
  chars = collect(line)
  n = length(chars)
  out = String[]
  pos = 1
  for (_, w) in fmt
    pos > n && break
    last = min(pos + w - 1, n)
    push!(out, strip(String(chars[pos:last])))
    pos += w
  end
  out
end

parse_int_field(s) = isempty(s) ? 0 : parse(Int, s)

function parse_real_field(s)
  isempty(s) && return 0.0
  t = replace(s, 'D' => 'E', 'd' => 'e')
  v = tryparse(Float64, t)
  if v === nothing
    # Fortran may drop the 'E' for 3-digit exponents: 1.234-100
    m = match(r"^([+-]?[0-9.]+)([+-]\d+)$", t)
    m === nothing && error("Cannot parse real number '$s'")
    v = parse(Float64, m.captures[1] * "E" * m.captures[2])
  end
  v
end

# Integers of a line written with an "(Ni W)" format. Falls back to
# whitespace splitting when the line does not look fixed-width.
function read_ints(line::AbstractString, fmt)
  vals = fixed_width_split(line, fmt)
  out = Int[]
  for v in vals
    isempty(v) && continue
    x = tryparse(Int, v)
    if x === nothing
      # not fixed width after all
      return [parse(Int, t) for t in split(line)]
    end
    push!(out, x)
  end
  out
end

header_fields(line) = [strip(f) for f in split(split(line, '!')[1], ',')]

# ---------------------------------------------------------------------------
# Main reader
# ---------------------------------------------------------------------------

"""
    read_cdb(filename) -> CDBMesh

Read nodes, elements, element types and components from an ANSYS `.cdb` file.
"""
function read_cdb(filename::AbstractString)
  lines = readlines(filename)
  read_cdb_lines(lines)
end

function read_cdb_lines(lines::Vector{String})
  node_ids = Int[]
  xyz = Float64[]
  elements = CDBElement[]
  etypes = Dict{Int,Int}()
  ncomps = Pair{String,Vector{Int}}[]
  ecomps = Pair{String,Vector{Int}}[]

  params = Dict{String,Int}()          # integer APDL parameters (*set,tid,4)
  esel = nothing                       # current element selection (nothing = all)
  body = nothing                       # Workbench body name of the next EBLOCK
  bcomps = Pair{String,Vector{Int}}[]  # one element component per Workbench body

  nl = length(lines)
  i = 1
  while i <= nl
    raw = rstrip(lines[i], '\r')
    up = uppercase(lstrip(raw))
    if startswith(up, "ET,")
      f = header_fields(raw)
      if length(f) >= 3
        id = apdl_int(f[2], params)
        m = match(r"(\d+)\s*$", f[3])
        if id !== nothing && m !== nothing
          etypes[id] = parse(Int, m.captures[1])
        end
      end
      i += 1
    elseif startswith(up, "*SET,")
      f = header_fields(raw)
      if length(f) >= 3
        v = tryparse(Int, f[3])
        v === nothing || (params[uppercase(f[2])] = v)
      end
      i += 1
    elseif (m = match(r"^\s*([A-Za-z_]\w*)\s*=\s*([+-]?\d+)\s*(!.*)?$", raw)) !== nothing
      params[uppercase(m.captures[1])] = parse(Int, m.captures[2])
      i += 1
    elseif startswith(up, "ESEL,")
      esel = apply_esel(esel, header_fields(raw), elements, etypes, params)
      i += 1
    elseif startswith(up, "CM,")
      f = header_fields(raw)
      if length(f) >= 3 && startswith(uppercase(f[3]), "ELEM")
        ids = esel === nothing ? [e.id for e in elements] : sort!(collect(esel))
        push!(ecomps, f[2] => ids)
      end
      i += 1
    elseif startswith(up, "NBLOCK")
      i = read_nblock!(node_ids, xyz, lines, i)
    elseif (m = match(r"Elements for Body\s+\d+\s+\"(.*)\"", raw)) !== nothing
      body = replace(m.captures[1], '\\' => '/')   # Workbench body name
      i += 1
    elseif startswith(up, "/WB,ELEM,END") || startswith(up, "/WB,CONTACT")
      body = nothing
      i += 1
    elseif startswith(up, "EBLOCK")
      n0 = length(elements)
      i = read_eblock!(elements, lines, i)
      if body !== nothing && length(elements) > n0
        k = findfirst(c -> first(c) == body, bcomps)
        k === nothing && (push!(bcomps, body => Int[]); k = length(bcomps))
        append!(last(bcomps[k]), (e.id for e in view(elements, n0+1:length(elements))))
      end
    elseif startswith(up, "CMBLOCK")
      i = read_cmblock!(ncomps, ecomps, lines, i)
    else
      i += 1
    end
  end

  # Workbench bodies go first; a named selection with the same name wins
  for b in reverse(bcomps)
    any(c -> first(c) == first(b), ecomps) || pushfirst!(ecomps, b)
  end
  coords = reshape(xyz, 3, :)
  CDBMesh(node_ids, coords, elements, etypes, ncomps, ecomps)
end

function read_nblock!(node_ids, xyz, lines, i)
  fmt = parse_fortran_format(lines[i+1])
  ni = count(f -> f[1] == :i, fmt)
  i += 2
  while i <= length(lines)
    line = rstrip(lines[i], '\r')
    s = strip(line)
    i += 1
    isempty(s) && continue
    startswith(uppercase(s), "N,") && break
    vals = fixed_width_split(line, fmt)
    id = parse_int_field(vals[1])
    id == -1 && break
    push!(node_ids, id)
    for k in 1:3
      j = ni + k
      push!(xyz, j <= length(vals) ? parse_real_field(vals[j]) : 0.0)
    end
  end
  i
end

function read_eblock!(elements, lines, i)
  f = header_fields(lines[i])
  solid = length(f) >= 3 && uppercase(f[3]) == "SOLID"
  fmt = parse_fortran_format(lines[i+1])
  i += 2
  while i <= length(lines)
    line = rstrip(lines[i], '\r')
    i += 1
    isempty(strip(line)) && continue
    v = read_ints(line, fmt)
    isempty(v) && continue
    v[1] == -1 && break
    if solid
      # mat type real sec esys birth solidref shape nnodes unused elnum nodes...
      length(v) < 11 && error("Malformed SOLID EBLOCK line: '$line'")
      nn = v[9]
      nodes = v[12:end]
      while length(nodes) < nn && i <= length(lines)
        push!(nodes, read_ints(rstrip(lines[i], '\r'), fmt)...)
        i += 1
      end
      resize!(nodes, nn)
      push!(elements, CDBElement(v[11], v[2], v[1], v[3], v[4], nodes))
    else
      # elnum type real mat sec nodes... (surface / contact elements)
      length(v) < 6 && error("Malformed EBLOCK line: '$line'")
      push!(elements, CDBElement(v[1], v[2], v[4], v[3], v[5], filter(!=(0), v[6:end])))
    end
  end
  i
end

function read_cmblock!(ncomps, ecomps, lines, i)
  f = header_fields(lines[i])
  name = f[2]
  kind = uppercase(f[3])
  nitems = parse(Int, f[4])
  fmt = parse_fortran_format(lines[i+1])
  i += 2
  raw = Int[]
  while length(raw) < nitems && i <= length(lines)
    line = rstrip(lines[i], '\r')
    i += 1
    isempty(strip(line)) && continue
    append!(raw, read_ints(line, fmt))
  end
  ids = Int[]
  for v in raw
    if v < 0
      isempty(ids) && error("Component $name starts with a negative range")
      append!(ids, (ids[end]+1):(-v))
    else
      push!(ids, v)
    end
  end
  if kind == "NODE"
    push!(ncomps, name => ids)
  elseif kind == "ELEM" || kind == "ELEMENT"
    push!(ecomps, name => ids)
  end
  i
end

# ---------------------------------------------------------------------------
# Minimal APDL support for element selections, as written by Workbench:
#   esel,s,type,,2 / esel,a,mat,,1,3 / esel,all / cm,MATA,elem
# ---------------------------------------------------------------------------

function apdl_int(s::AbstractString, params)
  v = tryparse(Int, s)
  v === nothing ? get(params, uppercase(strip(s)), nothing) : v
end

function apply_esel(sel, f, elements, etypes, params)
  typ = length(f) >= 2 ? uppercase(f[2]) : "S"
  typ in ("ALL", "") && return nothing
  typ == "NONE" && return Set{Int}()
  allids = Set(e.id for e in elements)
  current = sel === nothing ? allids : sel
  typ == "INVE" && return setdiff(allids, current)
  item = length(f) >= 3 ? uppercase(f[3]) : ""
  vmin = length(f) >= 5 ? apdl_int(f[5], params) : nothing
  vmax = length(f) >= 6 && !isempty(f[6]) ? apdl_int(f[6], params) : vmin
  vinc = length(f) >= 7 && !isempty(f[7]) ? apdl_int(f[7], params) : 1
  if vmin === nothing || !(item in ("TYPE", "MAT", "REAL", "SEC", "SECN", "ELEM", "ENAME"))
    @warn "Unsupported selection command ignored: $(join(f, ','))"
    return sel
  end
  vmax === nothing && (vmax = vmin)
  vinc === nothing && (vinc = 1)
  vals = Set(vmin:vinc:vmax)
  field(e) = item == "TYPE" ? e.etype : item == "MAT" ? e.mat : item == "REAL" ? e.real :
             item == "ELEM" ? e.id : item == "ENAME" ? get(etypes, e.etype, 0) : e.secnum
  hit = Set(e.id for e in elements if field(e) in vals)
  typ == "S" ? hit :
  typ == "R" ? intersect(current, hit) :
  typ == "A" ? union(current, hit) :
  typ == "U" ? setdiff(current, hit) :
  (@warn("Unsupported selection type '$typ' ignored"); sel)
end
