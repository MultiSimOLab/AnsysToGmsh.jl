module AnsysToGmshGridapGmshExt

using AnsysToGmsh
using AnsysToGmsh: GmshMesh
using GridapGmsh
using GridapGmsh: gmsh
import GridapGmsh.Gridap

# Load a GmshMesh into the (in-memory) Gmsh API model
function load_into_gmsh!(mesh::GmshMesh)
  gmsh.model.add("cdb")
  ents = sort(mesh.entities, by=e -> (e.dim, e.tag))
  for e in ents
    gmsh.model.addDiscreteEntity(e.dim, e.tag)
  end
  # physical groups: name and member entities
  for (d, ptag, name) in mesh.physical_names
    members = [e.tag for e in ents if e.dim == d && ptag in e.physicals]
    gmsh.model.addPhysicalGroup(d, members, ptag)
    gmsh.model.setPhysicalName(d, ptag, name)
  end
  X = mesh.coords
  for e in ents
    isempty(e.nodes) && continue
    gmsh.model.mesh.addNodes(e.dim, e.tag, e.nodes, vec(X[:, e.nodes]))
  end
  tag = 0
  for e in ents
    isempty(e.elements) && continue
    n = length(e.elements)
    tags = collect(tag+1:tag+n); tag += n
    gmsh.model.mesh.addElements(e.dim, e.tag, [e.etype], [tags], [reduce(vcat, e.elements)])
  end
  nothing
end

function AnsysToGmsh.CDBDiscreteModel(cdbfile::AbstractString; has_affine_map=nothing,
                                      orient_if_simplex=nothing, kwargs...)
  mesh = convert_cdb(read_cdb(cdbfile); kwargs...)
  AnsysToGmsh.CDBDiscreteModel(mesh; has_affine_map, orient_if_simplex)
end

function AnsysToGmsh.CDBDiscreteModel(mesh::GmshMesh; has_affine_map=nothing, orient_if_simplex=nothing)
  gmsh.initialize()
  try
    gmsh.option.setNumber("General.Terminal", 0)
    load_into_gmsh!(mesh)
    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()
    return GmshDiscreteModel(gmsh; has_affine_map, orient_if_simplex)
  finally
    gmsh.isInitialized() == 1 && gmsh.finalize()
  end
end

# DiscreteModelFromFile("model.cdb") also works (default conversion options)
Gridap.DiscreteModelFromFile(filename::AbstractString, ::Val{:cdb}) = AnsysToGmsh.CDBDiscreteModel(filename)

end
