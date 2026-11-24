# AnsysToGmsh.jl

Convert ANSYS archive files (`.cdb`, from Mechanical APDL `CDWRITE` or
exported from Workbench) into Gmsh 4.1 meshes (`.msh`) that load directly
with `GmshDiscreteModel` in [Gridap](https://github.com/gridap/Gridap.jl) /
[HyperFEM](https://github.com/MultiSimOLab/HyperFEM.jl).

Named selections survive the trip: every ANSYS component becomes a Gmsh
physical group, and therefore a Gridap tag you can use for Dirichlet
conditions, boundary integrals or material regions.

```julia
using AnsysToGmsh
cdb2msh("model.cdb")                       # writes model.msh

using HyperFEM                             # or: using Gridap, GridapGmsh
model = GmshDiscreteModel("model.msh")
Ω  = Triangulation(model, tags = "MATA")   # body / element named selection
Γ  = BoundaryTriangulation(model, tags = "TOP")
V  = TestFESpace(model, reffe, dirichlet_tags = ["Y0"])
```

No dependencies beyond Julia itself.

## Installation

```julia
] add https://github.com/MultiSimOLab/AnsysToGmsh.jl
```

## Command line

```bash
julia --project=path/to/AnsysToGmsh.jl path/to/AnsysToGmsh.jl/bin/cdb2msh.jl model.cdb
julia --project=path/to/AnsysToGmsh.jl path/to/AnsysToGmsh.jl/bin/cdb2msh.jl modelos/*.cdb --split
```

Options: `--split`, `--linear`, `--no-boundary`, `--no-domain`,
`--no-materials`, `--boundary-only`, `--quiet` (see `--help`).

## What you get in the `.msh`

| Physical group | Dimension | Comes from |
|---|---|---|
| `domain` | cells | all converted elements |
| `material_<id>` | cells | ANSYS material number (`MAT`) |
| body name (`Solid`, `Solid1`, …) | cells | Workbench bodies (`/com, Elements for Body …`) |
| element named selections (`MATA`, …) | cells | `CMBLOCK,…,ELEM` and Workbench `esel`/`cm,…,elem` |
| `boundary` | faces, edges, vertices | whole boundary of the mesh |
| node named selections (`TOP`, `Y0`, …) | faces, edges, vertices | `CMBLOCK,…,NODE` |

A face (or edge) belongs to a node named selection when all its corner nodes
do. Groups are written on every lower dimension (vertices, edges, faces), as
GridapGmsh expects, so a vertex on the edge shared by `TOP` and `Y0` is in
both tags. Interior faces are tagged too (e.g. an `INTERFACE` between two
bodies) unless `interior_faces = false`.

Node and element numbers are renumbered consecutively, unused nodes
(contact pilots, dropped midside nodes) are removed.

## Supported elements

| ANSYS | Gmsh / Gridap |
|---|---|
| SOLID185, 45, 5, 65, 69, 70, 96, 164, 278 | hex8 (degenerate forms: tet4, prism, pyramid) |
| SOLID186, 95, 90, 117, 122, 226, 231, 236, 279 | hex8 (midside nodes dropped); degenerate tet → tet10 |
| SOLID187, 92, 87, 98, 123, 227, 232, 237 | tet10 (tet4 with `order = 1`) |
| SOLID285, 72 | tet4 |
| PLANE182, 42, 55, 13 | quad4 / tri3 |
| PLANE183, 82, 77, 53, 121, 223, 230, 233 | quad4 (linearized) / tri6 |
| PLANE2, 35 | tri6 |

Contact, target, surface and shell elements (CONTA17x, TARGE170, SURF15x,
SHELL…) are skipped with a message.

Gridap only handles Lagrangian elements and a single element type per mesh,
which drives two choices:

* **Serendipity elements** (20-node bricks, 8-node quads) become linear.
  Quadratic Lagrangian elements (10-node tets, 6-node triangles) are kept.
  `order = 1` makes everything linear.
* **Mixed meshes** (Workbench *Hex Dominant*: bricks + prisms + pyramids +
  tets) raise an error by default. With `split = true` every cell is turned
  into linear tetrahedra in a conforming way: quadrilateral faces are cut
  along the diagonal through their smallest node id (so neighbours agree),
  pyramids into two tets, and bricks/prisms into one tet per face triangle
  around an interior point (the centroid, or the centre of the cell's kernel
  when the centroid does not give valid tets).

## API

```julia
cdb  = read_cdb("model.cdb")              # CDBMesh: raw nodes, elements, components
mesh = convert_cdb(cdb;
                   order = :auto,          # or 1
                   split = false,          # true / :always
                   boundary = "boundary",  # nothing to skip
                   domain = "domain",      # nothing to skip
                   materials = "material_",# nothing to skip
                   interior_faces = true,
                   verbose = true)
write_msh("model.msh", mesh)
```

## Exporting from Workbench

In Mechanical: *Model* (right click) → *Export* → *Ansys Input File* (or add a
command snippet with `cdwrite,db,model,cdb`). Named selections of bodies and
faces are written as components automatically. Use *Shared Topology* in the
geometry when bodies must share nodes (no contact).

## Tests

The test suite uses real Workbench exports (tet linear/quadratic, Hex
Dominant, a body with a hole and a non-contiguous named selection, and a model
with contact elements) and checks volumes, areas, tags and conformity, plus a
full round trip through `GmshDiscreteModel`.

```julia
] test AnsysToGmsh
```
