#!/usr/bin/env julia
# Command line converter.
#
#   julia --project=path/to/AnsysToGmsh.jl bin/cdb2msh.jl model.cdb [model.msh] [options]
#   julia --project=path/to/AnsysToGmsh.jl bin/cdb2msh.jl "modelos ANSYS/"*.cdb --split
#
# Options
#   --split            split mixed meshes (Hex Dominant) into linear tetrahedra
#   --linear           make every element linear (order = 1)
#   --no-boundary      do not add the automatic "boundary" group
#   --no-domain        do not add the automatic "domain" group
#   --no-materials     do not add "material_<id>" groups
#   --boundary-only    tag only boundary faces (ignore interior interfaces)
#   --quiet            no summary
using AnsysToGmsh

function main(args)
  opts = filter(startswith("--"), args)
  files = filter(!startswith("--"), args)
  if isempty(files) || "--help" in opts
    println(read(@__FILE__, String)[findfirst("# Command", read(@__FILE__, String))[1]:findfirst("using AnsysToGmsh", read(@__FILE__, String))[1]-1])
    return
  end
  kw = (; split="--split" in opts,
          order="--linear" in opts ? 1 : :auto,
          boundary="--no-boundary" in opts ? nothing : "boundary",
          domain="--no-domain" in opts ? nothing : "domain",
          materials="--no-materials" in opts ? nothing : "material_",
          interior_faces=!("--boundary-only" in opts),
          verbose=!("--quiet" in opts))
  pairs = if length(files) == 2 && endswith(lowercase(files[2]), ".msh")
    [(files[1], files[2])]
  else
    [(f, splitext(f)[1] * ".msh") for f in files]
  end
  for (cdb, msh) in pairs
    println("→ ", cdb)
    try
      cdb2msh(cdb, msh; kw...)
      println("  wrote ", msh)
    catch e
      println("  FAILED: ", sprint(showerror, e))
    end
  end
end

main(ARGS)
