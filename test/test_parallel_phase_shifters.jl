# T68: a phase shifter sharing its arc with another branch. The direct map holds one branch per
# arc, so before, the arc's later branch overwrote the earlier one (a Line beside a phase shifter,
# or two phase shifters): the Y-bus was right but the overwritten branch had no forward entry.
# Now the arc's branches form a parallel set, each keeping its own entry.
function _pps_system(; α12 = 0.1, α23 = (0.05, -0.08))
    sys = System(100.0)
    buses = ACBus[]
    for i in 1:3
        b = ACBus(; number = i, name = "pps$i", available = true,
            bustype = i == 1 ? ACBusTypes.REF : ACBusTypes.PQ, angle = 0.0,
            magnitude = 1.0, voltage_limits = (min = 0.9, max = 1.1), base_voltage = 230.0)
        add_component!(sys, b)
        push!(buses, b)
    end
    pst(name, f, t, α, x) = PhaseShiftingTransformer(; name = name, available = true,
        active_power_flow = 0.0, reactive_power_flow = 0.0,
        arc = Arc(; from = buses[f], to = buses[t]), r = 0.01, x = x,
        primary_shunt = 0.0, tap = 1.0, α = α, rating = 1.0, base_power = 100.0)
    line(name, f, t, x) = Line(; name = name, available = true, active_power_flow = 0.0,
        reactive_power_flow = 0.0, arc = Arc(; from = buses[f], to = buses[t]),
        r = 0.01, x = x, b = (from = 0.01, to = 0.01), rating = 1.0,
        angle_limits = (min = -π / 2, max = π / 2))
    for br in (line("L12", 1, 2, 0.1), pst("P12", 1, 2, α12, 0.2),
        pst("P23a", 2, 3, α23[1], 0.15), pst("P23b", 2, 3, α23[2], 0.25),
        line("L13", 1, 3, 0.3))
        add_component!(sys, br)
    end
    return sys
end

@testset "T68: a phase shifter in parallel keeps its own branch entry" begin
    sys = _pps_system()
    ybus = Ybus(sys)
    nr = ybus.network_reduction_data
    direct = PNM.get_direct_branch_map(nr)
    parallel = PNM.get_parallel_branch_map(nr)
    # A line beside a phase shifter (two types), and two phase shifters.
    @test Set(PSY.get_name.(parallel[(1, 2)].branches)) == Set(["L12", "P12"])
    @test parallel[(1, 2)] isa PNM.MixedBranchesParallel
    @test Set(PSY.get_name.(parallel[(2, 3)].branches)) == Set(["P23a", "P23b"])
    @test !haskey(direct, (1, 2)) && !haskey(direct, (2, 3))
    @test PSY.get_name(direct[(1, 3)]) == "L13"
    # Every branch is in exactly one forward map and its reverse map.
    rev_d = PNM.get_reverse_direct_branch_map(nr)
    rev_p = PNM.get_reverse_parallel_branch_map(nr)
    for br in PSY.get_components(PSY.ACTransmission, sys)
        @test haskey(rev_d, br) ⊻ haskey(rev_p, br)
    end
    @test length(rev_d) == length(direct) == 1
    @test length(rev_p) == 4
    # The Y-bus is the sum of the branches' own blocks, whatever the maps.
    bl = PNM.get_bus_lookup(ybus)
    Y = zeros(ComplexF64, 3, 3)
    for br in PSY.get_components(PSY.ACTransmission, sys)
        y11, y12, y21, y22 = PNM.ybus_branch_entries(br, nr)
        f, t = PNM.get_arc_tuple(br)
        i, j = bl[f], bl[t]
        Y[i, i] += y11
        Y[i, j] += y12
        Y[j, i] += y21
        Y[j, j] += y22
    end
    @test Matrix(ybus.data) ≈ Y atol = 1e-12
    # Typed maps list the mixed set under both types.
    PNM.populate_branch_maps_by_type!(nr)
    by_type = nr.all_branch_maps_by_type.parallel_branch_map
    @test haskey(by_type[PSY.Line], (1, 2))
    @test haskey(by_type[PSY.PhaseShiftingTransformer], (1, 2))
    @test haskey(by_type[PSY.PhaseShiftingTransformer], (2, 3))
    @test nr.component_to_reduction_name_map[PSY.PhaseShiftingTransformer]["P12"] ==
          PNM.get_name(parallel[(1, 2)])
    # A physical equivalent exists only where the phase shifts agree.
    @test_throws ErrorException PNM.get_equivalent_physical_branch_parameters(
        parallel[(2, 3)], nr)
    eq_nr = Ybus(_pps_system(; α23 = (0.05, 0.05))).network_reduction_data
    eq = PNM.get_equivalent_physical_branch_parameters(
        PNM.get_parallel_branch_map(eq_nr)[(2, 3)], eq_nr)
    @test abs(PNM.get_equivalent_shift(eq)) ≈ 0.05 atol = 1e-9
end
