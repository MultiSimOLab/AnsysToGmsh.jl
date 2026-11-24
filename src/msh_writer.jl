# Writer for Gmsh MSH 4.1 (ASCII), the format read by GridapGmsh.

"""
    write_msh(filename, mesh::GmshMesh)
    write_msh(io::IO, mesh::GmshMesh)

Write `mesh` as an ASCII Gmsh 4.1 file.
"""
write_msh(filename::AbstractString, mesh::GmshMesh) = open(io -> write_msh(io, mesh), filename, "w")

function write_msh(io::IO, mesh::GmshMesh)
  X = mesh.coords
  ents = sort(mesh.entities, by=e -> (e.dim, e.tag))

  println(io, "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat")

  if !isempty(mesh.physical_names)
    println(io, "\$PhysicalNames")
    println(io, length(mesh.physical_names))
    for (d, t, name) in sort(mesh.physical_names)
      println(io, d, " ", t, " \"", name, "\"")
    end
    println(io, "\$EndPhysicalNames")
  end

  # Entities
  println(io, "\$Entities")
  println(io, join([count(e -> e.dim == d, ents) for d in 0:3], " "))
  for e in ents
    phys = string(length(e.physicals), (" $p" for p in e.physicals)...)
    if e.dim == 0
      n = isempty(e.nodes) ? e.elements[1][1] : e.nodes[1]
      println(io, e.tag, " ", fmt(X[1, n]), " ", fmt(X[2, n]), " ", fmt(X[3, n]), " ", phys)
    else
      lo, hi = bbox(X, e)
      println(io, e.tag, " ", join(fmt.(lo), " "), " ", join(fmt.(hi), " "), " ", phys, " 0")
    end
  end
  println(io, "\$EndEntities")

  # Nodes, one block per entity that owns nodes
  nblocks = count(e -> !isempty(e.nodes), ents)
  nn = size(X, 2)
  println(io, "\$Nodes")
  println(io, nblocks, " ", nn, " 1 ", nn)
  for e in ents
    isempty(e.nodes) && continue
    println(io, e.dim, " ", e.tag, " 0 ", length(e.nodes))
    for n in e.nodes
      println(io, n)
    end
    for n in e.nodes
      println(io, fmt(X[1, n]), " ", fmt(X[2, n]), " ", fmt(X[3, n]))
    end
  end
  println(io, "\$EndNodes")

  # Elements, consecutive tags dimension by dimension
  eblocks = count(e -> !isempty(e.elements), ents)
  ne = sum(e -> length(e.elements), ents; init=0)
  println(io, "\$Elements")
  println(io, eblocks, " ", ne, " 1 ", ne)
  tag = 0
  for e in ents
    isempty(e.elements) && continue
    println(io, e.dim, " ", e.tag, " ", e.etype, " ", length(e.elements))
    for el in e.elements
      tag += 1
      print(io, tag)
      for n in el
        print(io, " ", n)
      end
      println(io)
    end
  end
  println(io, "\$EndElements")
  nothing
end

fmt(x::Float64) = string(x)

function bbox(X, e::MshEntity)
  lo = fill(Inf, 3)
  hi = fill(-Inf, 3)
  for el in e.elements, n in el, k in 1:3
    lo[k] = min(lo[k], X[k, n])
    hi[k] = max(hi[k], X[k, n])
  end
  lo, hi
end
