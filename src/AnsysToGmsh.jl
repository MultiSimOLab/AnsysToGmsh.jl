"""
    AnsysToGmsh

Convert ANSYS archive files (`.cdb`) into Gmsh 4.1 meshes (`.msh`) that can
be loaded with `GmshDiscreteModel` (GridapGmsh) and used in HyperFEM.

Main entry point: [`cdb2msh`](@ref). Lower level: [`read_cdb`](@ref),
[`convert_cdb`](@ref), [`write_msh`](@ref).
"""
module AnsysToGmsh

export read_cdb, convert_cdb, write_msh, cdb2msh, supported_ansys_elements
export CDBMesh, CDBElement, GmshMesh, MshEntity

include("cdb_reader.jl")
include("elements.jl")
include("conversion.jl")
include("msh_writer.jl")

"""
    cdb2msh(cdbfile, mshfile=splitext(cdbfile)[1] * ".msh"; kwargs...) -> mshfile

Read an ANSYS `.cdb` file and write a Gmsh 4.1 `.msh` file ready for
`GmshDiscreteModel`. Keyword arguments are those of [`convert_cdb`](@ref).

```julia
using AnsysToGmsh
cdb2msh("model.cdb")                 # writes model.msh

using HyperFEM
model = GmshDiscreteModel("model.msh")
```
"""
function cdb2msh(cdbfile::AbstractString, mshfile::AbstractString=splitext(cdbfile)[1] * ".msh"; kwargs...)
  cdb = read_cdb(cdbfile)
  mesh = convert_cdb(cdb; kwargs...)
  write_msh(mshfile, mesh)
  mshfile
end

end # module
