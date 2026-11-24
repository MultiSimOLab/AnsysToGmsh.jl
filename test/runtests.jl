using AnsysToGmsh
using Test

const A = AnsysToGmsh
const DATA = joinpath(@__DIR__, "data")
const TMP = mktempdir()

# ---------------------------------------------------------------------------
# helpers: measures computed from the GmshMesh itself (linear geometry)
# ---------------------------------------------------------------------------
pt(m, n) = m.coords[:, n]
cross3(u, v) = [u[2]v[3] - u[3]v[2], u[3]v[1] - u[1]v[3], u[1]v[2] - u[2]v[1]]
tetvol(m, e) = abs(A.orient3(pt(m, e[1]), pt(m, e[2]), pt(m, e[3]), pt(m, e[4]))) / 6
trisurf(m, e) = sqrt(sum(abs2, cross3(pt(m, e[2]) - pt(m, e[1]), pt(m, e[3]) - pt(m, e[1])))) / 2
hexvol(m, e) = sum(tetvol(m, [e[i] for i in t]) for t in ((1, 2, 4, 5), (2, 3, 4, 7), (2, 5, 6, 7), (4, 5, 7, 8), (2, 4, 5, 7)))
measure(m, e, t) = t in (A.GMSH_TET4, A.GMSH_TET10) ? tetvol(m, e) :
                   t == A.GMSH_HEX8 ? hexvol(m, e) :
                   t in (A.GMSH_TRI3, A.GMSH_TRI6) ? trisurf(m, e) :
                   t == A.GMSH_QUAD4 ? trisurf(m, e[[1, 2, 3]]) + trisurf(m, e[[1, 3, 4]]) : error()

function group_measure(m, dim, name)
  tags = [t for (d, t, n) in m.physical_names if d == dim && n == name]
  @assert length(tags) == 1 "group $name (dim $dim) not found"
  tot = 0.0
  for e in m.entities
    (e.dim == dim && tags[1] in e.physicals) || continue
    tot += sum(el -> measure(m, el, e.etype), e.elements)
  end
  tot
end

group_nodes(m, dim, name) = begin
  tags = [t for (d, t, n) in m.physical_names if d == dim && n == name]
  Set(n for e in m.entities if e.dim == dim && tags[1] in e.physicals for el in e.elements for n in el)
end

@testset "AnsysToGmsh" begin

  @testset "Fortran formats" begin
    @test A.parse_fortran_format("(3i9,6e21.13e3)") == [fill((:i, 9), 3); fill((:r, 21), 6)]
    @test A.parse_fortran_format("(1i9,3e20.9e3)") == [(:i, 9); fill((:r, 20), 3)]
    @test A.parse_real_field("1.000000000E+02") == 100.0
    @test A.parse_real_field("1.5-100") == 1.5e-100
  end

  @testset "element conversion" begin
    # SOLID187: ANSYS I J K L M N O P Q R -> Gmsh tet10
    r = A.ansys_to_gmsh(187, collect(1:10))
    @test r[1] == A.GMSH_TET10 && r[2] == [1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
    @test A.ansys_to_gmsh(187, collect(1:10); order=1)[1:2] == (A.GMSH_TET4, [1, 2, 3, 4])
    # SOLID185 brick and its degenerate forms
    @test A.ansys_to_gmsh(185, collect(1:8))[1:2] == (A.GMSH_HEX8, collect(1:8))
    @test A.ansys_to_gmsh(185, [1, 2, 3, 3, 4, 4, 4, 4])[1:2] == (A.GMSH_TET4, [1, 2, 3, 4])
    @test A.ansys_to_gmsh(185, [1, 2, 3, 4, 5, 5, 5, 5])[1:2] == (A.GMSH_PYR5, [1, 2, 3, 4, 5])
    @test A.ansys_to_gmsh(185, [1, 2, 3, 3, 4, 5, 6, 6])[1:2] == (A.GMSH_PRISM6, [1, 2, 3, 4, 5, 6])
    # SOLID186: bricks are linearized, degenerate tets stay quadratic
    r = A.ansys_to_gmsh(186, collect(1:20))
    @test r == (A.GMSH_HEX8, collect(1:8), :linearized)
    n = [1, 2, 3, 3, 4, 4, 4, 4, 11, 12, 3, 14, 4, 4, 4, 4, 17, 18, 19, 19]
    @test A.ansys_to_gmsh(186, n)[1:2] == (A.GMSH_TET10, [1, 2, 3, 4, 11, 12, 14, 17, 19, 18])
    # 2-D
    @test A.ansys_to_gmsh(182, [1, 2, 3, 4])[1:2] == (A.GMSH_QUAD4, [1, 2, 3, 4])
    @test A.ansys_to_gmsh(182, [1, 2, 3, 3])[1:2] == (A.GMSH_TRI3, [1, 2, 3])
    @test A.ansys_to_gmsh(183, [1, 2, 3, 3, 5, 6, 3, 8])[1:2] == (A.GMSH_TRI6, [1, 2, 3, 5, 6, 8])
    @test A.ansys_to_gmsh(174, collect(1:8)) === nothing
  end

  @testset "small hand-written CDB (2 bricks, ranges, esel/cm)" begin
    cdb = """
    et,1,185
    nblock,3,,12
    (1i9,3e20.9e3)
            1     0.000000000E+00     0.000000000E+00     0.000000000E+00
            2     1.000000000E+00     0.000000000E+00     0.000000000E+00
            3     1.000000000E+00     1.000000000E+00     0.000000000E+00
            4     0.000000000E+00     1.000000000E+00     0.000000000E+00
            5     0.000000000E+00     0.000000000E+00     1.000000000E+00
            6     1.000000000E+00     0.000000000E+00     1.000000000E+00
            7     1.000000000E+00     1.000000000E+00     1.000000000E+00
            8     0.000000000E+00     1.000000000E+00     1.000000000E+00
            9     2.000000000E+00     0.000000000E+00     0.000000000E+00
           10     2.000000000E+00     1.000000000E+00     0.000000000E+00
           11     2.000000000E+00     0.000000000E+00     1.000000000E+00
           12     2.000000000E+00     1.000000000E+00     1.000000000E+00
    N,R5.3,LOC,       -1,
    eblock,19,solid,,2
    (19i9)
            1        1        1        1        0        0        0        0        8        0        1        1        2        3        4        5        6        7        8
            2        1        1        1        0        0        0        0        8        0        2        2        9       10        3        6       11       12        7
           -1
    esel,s,mat,,2
    cm,RIGHT,elem
    esel,all
    CMBLOCK,XMAX    ,NODE,        2
    (8i10)
             9       -12
    CMBLOCK,X0      ,NODE,        4
    (8i10)
             1         4         5         8
    """
    f = joinpath(TMP, "two_bricks.cdb")
    write(f, cdb)
    c = read_cdb(f)
    @test length(c.node_ids) == 12 && length(c.elements) == 2
    @test c.node_components[1] == ("XMAX" => [9, 10, 11, 12])
    @test c.element_components[1] == ("RIGHT" => [2])
    m = convert_cdb(c; verbose=false)
    @test group_measure(m, 3, "domain") ≈ 2
    @test group_measure(m, 3, "material_2") ≈ 1
    @test group_measure(m, 3, "RIGHT") ≈ 1
    @test group_measure(m, 2, "XMAX") ≈ 1
    @test group_measure(m, 2, "X0") ≈ 1
    @test group_measure(m, 2, "boundary") ≈ 10
    @test size(m.coords, 2) == 12
    out = cdb2msh(f; verbose=false)
    @test isfile(out) && startswith(read(out, String), "\$MeshFormat\n4.1 0 8")
  end

  @testset "Workbench tet meshes (SOLID185 degenerate, SOLID187)" begin
    for (file, T) in (("model2_tet_order1.cdb", A.GMSH_TET4), ("model2_tet_order2.cdb", A.GMSH_TET10))
      c = read_cdb(joinpath(DATA, file))
      @test first.(c.element_components) == ["Solid", "Solid1", "MATA", "MATB", "ALL_SOLID"]
      @test first.(c.node_components) == ["Y0", "Y1", "INTERFACE", "TOP", "EDGE_Y", "ORIGIN"]
      m = convert_cdb(c; verbose=false)
      @test all(e -> e.dim != 3 || e.etype == T, m.entities)
      @test group_measure(m, 3, "ALL_SOLID") ≈ 1e6
      @test group_measure(m, 3, "MATA") ≈ 5e5
      @test group_measure(m, 3, "MATB") ≈ 5e5
      for g in ("Y0", "Y1", "TOP", "INTERFACE")
        @test group_measure(m, 2, g) ≈ 1e4
      end
      @test group_measure(m, 2, "boundary") ≈ 6e4
      # a vertex on the edge shared by TOP and Y0 belongs to both groups
      @test !isempty(intersect(group_nodes(m, 0, "TOP"), group_nodes(m, 0, "Y0")))
      @test length(group_nodes(m, 0, "ORIGIN")) == 1
      # node tags consecutive in file order
      @test sort(vcat([e.nodes for e in m.entities]...)) == 1:size(m.coords, 2)
    end
  end

  @testset "non-adjacent faces in one named selection" begin
    m = convert_cdb(read_cdb(joinpath(DATA, "model5_notadjacement_tet_order2.cdb")); verbose=false)
    @test group_measure(m, 2, "NOTADJACEMENT") ≈ 1e4 rtol = 1e-8
    @test group_measure(m, 3, "ALL_SOLID") ≈ group_measure(m, 3, "MATA") + group_measure(m, 3, "MATB")
  end

  @testset "contact elements are ignored" begin
    c = read_cdb(joinpath(DATA, "model_contact.cdb"))
    @test sort(unique(values(c.etypes))) == [170, 174, 185]
    m = convert_cdb(c; verbose=false)
    @test group_measure(m, 3, "ALL_SOLID") ≈ 1e6
    @test group_measure(m, 2, "boundary") ≈ 8e4      # both sides of the gap
  end

  @testset "Hex Dominant: error, then conforming split" begin
    c = read_cdb(joinpath(DATA, "model3_hexdominant_order2.cdb"))
    @test_throws ErrorException convert_cdb(c; verbose=false)
    m = convert_cdb(c; split=true, verbose=false)
    @test all(e -> e.dim != 3 || e.etype == A.GMSH_TET4, m.entities)
    for e in m.entities, el in (e.dim == 3 ? e.elements : ())
      @test A.orient3((pt(m, n) for n in el)...) > 0
    end
    @test group_measure(m, 3, "domain") ≈ 1e6
    @test group_measure(m, 3, "Solid") ≈ 5e5
    @test group_measure(m, 2, "boundary") ≈ 6e4
    @test group_measure(m, 2, "INTERFACE") ≈ 1e4
    # conformity: every interior triangle shared by exactly two tets
    faces = Dict{Vector{Int},Int}()
    for e in m.entities, el in (e.dim == 3 ? e.elements : ())
      for f in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4))
        k = sort(el[collect(f)]); faces[k] = get(faces, k, 0) + 1
      end
    end
    @test maximum(values(faces)) == 2
    @test sum(k -> trisurf(m, k), (k for (k, v) in faces if v == 1)) ≈ 6e4
  end

  # Full round trip through GridapGmsh, as HyperFEM does
  if Base.find_package("GridapGmsh") !== nothing
    @eval using Gridap, GridapGmsh, Gridap.Geometry
    @testset "GridapGmsh round trip" begin
      for (file, kw) in (("model2_tet_order2.cdb", (;)), ("model3_hexdominant_order2.cdb", (; split=true)))
        msh = cdb2msh(joinpath(DATA, file), joinpath(TMP, replace(file, ".cdb" => ".msh")); verbose=false, kw...)
        model = GmshDiscreteModel(msh)
        vol(t) = sum(∫(1)Measure(Triangulation(model, tags=t), 4))
        area(t) = sum(∫(1)Measure(BoundaryTriangulation(model, tags=t), 4))
        @test vol("MATA") ≈ 5e5
        @test vol("ALL_SOLID") ≈ 1e6
        @test area("Y0") ≈ 1e4
        @test area("TOP") ≈ 1e4
        @test area("boundary") ≈ 6e4
        lab = get_face_labeling(model)
        @test any(get_face_mask(lab, "TOP", 0) .& get_face_mask(lab, "Y0", 0))
        # in-memory route gives the very same model
        mem = CDBDiscreteModel(joinpath(DATA, file); verbose=false, kw...)
        lab2 = get_face_labeling(mem)
        @test get_tag_name(lab2) == get_tag_name(lab)
        @test get_cell_node_ids(mem) == get_cell_node_ids(model)
        @test Gridap.Geometry.get_node_coordinates(get_grid(mem)) == Gridap.Geometry.get_node_coordinates(get_grid(model))
        for t in ("MATA", "TOP", "Y0", "INTERFACE", "boundary"), d in 0:3
          @test get_face_mask(lab2, t, d) == get_face_mask(lab, t, d)
        end
      end
      @test num_cells(DiscreteModelFromFile(joinpath(DATA, "model2_tet_order1.cdb"))) == 1619
    end
  end
end
