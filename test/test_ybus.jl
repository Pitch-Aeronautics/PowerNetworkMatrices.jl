@testset "Test Ybus Matrix" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    buses_5 = nodes5()
    branches_5 = branches5(buses_5)
    sys_5 = System(100.0)
    for b in buses_5
        add_component!(sys_5, b)
    end
    for br in branches_5
        add_component!(sys_5, br)
    end
    buses_14 = nodes14()
    branches_14 = branches14(buses_14)
    sys_14 = System(100.0)
    for b in buses_14
        add_component!(sys_14, b)
    end
    for br in branches_14
        add_component!(sys_14, br)
    end
    Ybus5 = Ybus(sys_5)

    I, J, V = findnz(Ybus5.data)
    indices = collect(zip(I, J))

    for i in indices
        @test isapprox(Ybus5[i[1], i[2]], Ybus5_matpower[i[1], i[2]], atol = 1e-2)
    end
    Ybus14 = Ybus(sys_14)
    I, J, V = findnz(Ybus14.data)
    indices = collect(zip(I, J))

    for i in indices
        @test isapprox(Ybus14[i[1], i[2]], Ybus14_matpower[i[1], i[2]], atol = 1e-2)
    end

    bus5 = get_component(Bus, sys, "nodeE")
    set_number!(bus5, 10)
    Y5NS = Ybus(sys)
    @test isapprox(Y5NS[10, 4], -3.3336667 + 33.336667im, atol = 1e-4)

    Y5NS = Ybus(sys_5)
    for buf in Bu5NS_ids
        for but in Bu5NS_ids
            @test isapprox(Y5NS[buf, but], Ybus5_matpower[buf, but], atol = 1e-3)
        end
    end

    @test isapprox(
        Ybus5[buses_5[1], buses_5[2]],
        (-3.5234840209999647 + 35.234840209999646im),
        atol = 1e-4,
    )

    t_sys5_re = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    # Make 2 islands. Island 1: 1 - 5. Island 2: 2, 3 ,4
    remove_component!(Line, t_sys5_re, "1")
    remove_component!(Line, t_sys5_re, "2")
    remove_component!(Line, t_sys5_re, "6")
    @test length(Ybus(t_sys5_re).subnetwork_axes) == 2

    t2_sys5_re = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    # Remove lines. Don't cause islands
    remove_component!(Line, t2_sys5_re, "3")
    remove_component!(Line, t2_sys5_re, "5")
    @test length(Ybus(t2_sys5_re).subnetwork_axes) == 1

    sys_3bus = PSB.build_system(PSB.PSYTestSystems, "psse_3bus_gen_cls_sys")
    bus_103 = PSY.get_component(PSY.ACBus, sys_3bus, "BUS 3")
    fix_shunt = FixedAdmittance("FixAdm_Bus3", true, bus_103, 0.0 + 0.2im)
    add_component!(sys_3bus, fix_shunt)
    Ybus3 = Ybus(sys_3bus)
    I, J, V = findnz(Ybus3.data)
    indices = collect(zip(I, J))
    for i in indices
        @test isapprox(Ybus3.data[i[1], i[2]], Ybus3_matpower[i[1], i[2]], atol = 1e-5)
    end
end

@testset "Test Ybus with single bus not numbered 1" begin
    sys_single = System(100.0)
    bus = ACBus(;
        number = 5,
        name = "Bus5",
        available = true,
        bustype = ACBusTypes.REF,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = (min = 0.9, max = 1.05),
        base_voltage = 100.0,
    )
    add_component!(sys_single, bus)
    ybus_single = Ybus(sys_single)
    @test size(ybus_single) == (1, 1)
    @test PNM.get_ref_bus(ybus_single) == [5]
end

@testset "Test modification of units base when constructing Ybus" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    PSY.set_units_base_system!(sys, "NATURAL_UNITS")
    ybus = @test_logs (
        :warn,
        r"Setting the system unit base from NATURAL_UNITS to SYSTEM_BASE for matrix construction",
    ) match_mode = :any Ybus(sys)
end

@testset "Ybus stores ComplexF64: a near-zero-impedance branch is exact to 1e-9" begin
    # x = 1e-4 pu gives |y| ≈ 1e4 pu. In ComplexF32 the stored entry would carry a
    # relative error of ~6e-8 (≈ 6e-4 pu, 0.06 MVA at 1 pu voltage), larger than a
    # 0.01 MVA mismatch bar; in ComplexF64 it is exact to rounding. r > 0 keeps the
    # branch out of the always-on zero-impedance reduction (which needs r == 0).
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    line = get_component(Line, sys, "1")
    from_no = get_number(get_from(get_arc(line)))
    to_no = get_number(get_to(get_arc(line)))
    y_old = 1 / complex(get_r(line), get_x(line))
    yb_old = Ybus(sys)
    @test eltype(yb_old.data) == ComplexF64
    @test YBUS_ELTYPE == ComplexF64

    r, x = 1e-6, 1e-4
    set_r!(line, r)
    set_x!(line, x)
    y_new = 1 / complex(r, x)
    yb = Ybus(sys)
    @test eltype(yb.data) == ComplexF64
    @test size(yb) == size(yb_old)
    # Off-diagonal entries: −y for an untapped line with no parallel circuit.
    @test isapprox(yb[from_no, to_no], -y_new; rtol = 1e-9)
    @test isapprox(yb[to_no, from_no], -y_new; rtol = 1e-9)
    # Diagonal entries: the change is exactly the series admittance's change.
    for bus_no in (from_no, to_no)
        expected = yb_old[bus_no, bus_no] - y_old + y_new
        @test isapprox(yb[bus_no, bus_no], expected; rtol = 1e-9)
    end
    # The same entry rounded to ComplexF32 misses the bar, so the test discriminates.
    @test !isapprox(ComplexF64(ComplexF32(-y_new)), -y_new; rtol = 1e-9)
end
