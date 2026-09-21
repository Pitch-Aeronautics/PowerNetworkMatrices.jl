# ── Impedance correction (PSS/E transformer impedance correction tables) ─────
#
# Backport of Sienna-Platform/PowerNetworkMatrices.jl#298 (`psy6`, b49a158) onto v0.24.5.
# PowerSystems 5 parses a PSS/E table into an `ImpedanceCorrectionData` supplemental
# attribute on the transformer; `Ybus` scales the transformer's series impedance by the
# table's factor at the present tap ratio or phase-shift angle.

"""
Linearly interpolate `curve` at `x`, holding the end values outside the tabulated range.

PSS/E impedance correction tables are sparse and the operating point regularly sits outside
them; clamping matches PSS/E, which applies the nearest tabulated factor rather than
extrapolating.
"""
function _interpolate_correction_factor(curve::IS.PiecewiseLinearData, x::Real)
    points = IS.get_points(curve)
    x = clamp(x, points[1].x, points[end].x)
    for i in 1:(length(points) - 1)
        if x <= points[i + 1].x
            dx = points[i + 1].x - points[i].x
            iszero(dx) && return points[i].y
            t = (x - points[i].x) / dx
            return points[i].y + t * (points[i + 1].y - points[i].y)
        end
    end
    return points[end].y
end

# The operating point a table is read at: the tap ratio, or the phase shift in degrees (the
# tables' x-values are degrees; PowerSystems stores α in radians). On PowerSystems 5 a
# `Transformer2W` has no tap (ratio 1), and the α of a `Transformer2W`/`TapTransformer` is
# its winding-group shift.
_correction_tap(br::PSY.TwoWindingTransformer) = _get_tap(br)
_correction_angle(br::PSY.TwoWindingTransformer) = PSY.get_α(br)

function _correction_tap(w::ThreeWindingTransformerWinding)
    br = get_transformer(w)
    n = get_winding_number(w)
    n == 1 && return PSY.get_primary_turns_ratio(br)
    n == 2 && return PSY.get_secondary_turns_ratio(br)
    return PSY.get_tertiary_turns_ratio(br)
end

function _correction_angle(w::ThreeWindingTransformerWinding)
    br = get_transformer(w)
    n = get_winding_number(w)
    n == 1 && return PSY.get_α_primary(br)
    n == 2 && return PSY.get_α_secondary(br)
    return PSY.get_α_tertiary(br)
end

"""
The series-impedance multiplier `ict` prescribes for a transformer (or one winding of a
three-winding transformer) at its present operating point.
"""
function _evaluate_correction_table(
    br::Union{PSY.TwoWindingTransformer, ThreeWindingTransformerWinding},
    ict::PSY.ImpedanceCorrectionData,
)
    mode = PSY.get_transformer_control_mode(ict)
    if mode == PSY.ImpedanceCorrectionTransformerControlMode.TAP_RATIO
        x = abs(_correction_tap(br))
    else
        x = rad2deg(_correction_angle(br))
    end
    return _interpolate_correction_factor(PSY.get_impedance_correction_curve(ict), x)
end

# `WindingCategory` encodes the winding position directly (`TR2W_WINDING = 0`, then
# `PRIMARY_WINDING`/`SECONDARY_WINDING`/`TERTIARY_WINDING` = 1/2/3, matching
# `ThreeWindingTransformerWinding`'s winding number), so the code doubles as the map key.
_winding_index(category::PSY.WindingCategory) = Int(category.value)

"""
    build_impedance_correction_factors!(nr::NetworkReductionData, sys::PSY.System)

Evaluate every `PSY.ImpedanceCorrectionData` attached to a transformer in `sys` and cache the
resulting series-impedance multipliers on `nr`, keyed by `(transformer UUID, winding)`.

Evaluated once per `Ybus` build rather than per branch: resolving a component's supplemental
attributes queries the association store, which is far too costly for the assembly loop.
"""
function build_impedance_correction_factors!(nr::NetworkReductionData, sys::PSY.System)
    isempty(PSY.get_supplemental_attributes(PSY.ImpedanceCorrectionData, sys)) && return
    factors = nr.impedance_correction_factors
    # A 2W transformer has one circuit, so its table's winding tag carries no information and
    # is not validated by PSY; key on the 2W code regardless of the tag.
    for pair in PSY.get_component_supplemental_attribute_pairs(
        PSY.TwoWindingTransformer,
        PSY.ImpedanceCorrectionData,
        sys,
    )
        transformer = pair.component
        factors[(
            IS.get_uuid(transformer),
            _winding_index(PSY.WindingCategory.TR2W_WINDING),
        )] = _evaluate_correction_table(transformer, pair.supplemental_attribute)
    end
    for pair in PSY.get_component_supplemental_attribute_pairs(
        PSY.ThreeWindingTransformer,
        PSY.ImpedanceCorrectionData,
        sys,
    )
        transformer = pair.component
        ict = pair.supplemental_attribute
        winding = _winding_index(PSY.get_transformer_winding(ict))
        if !(winding in 1:3)
            error(
                "ImpedanceCorrectionData table $(PSY.get_table_number(ict)) on " *
                "ThreeWindingTransformer $(PSY.get_name(transformer)) is tagged " *
                "$(PSY.get_transformer_winding(ict)); expected PRIMARY_, SECONDARY_ or " *
                "TERTIARY_WINDING.",
            )
        end
        factors[(IS.get_uuid(transformer), winding)] = _evaluate_correction_table(
            ThreeWindingTransformerWinding(transformer, winding),
            ict,
        )
    end
    return
end

# Only transformers carry correction tables; every other branch kind is uncorrected.
_impedance_correction_factor(::PSY.ACTransmission, ::NetworkReductionData) = 1.0

_impedance_correction_factor(br::PSY.TwoWindingTransformer, nr::NetworkReductionData) =
    get(
        nr.impedance_correction_factors,
        (IS.get_uuid(br), _winding_index(PSY.WindingCategory.TR2W_WINDING)),
        1.0,
    )

_impedance_correction_factor(w::ThreeWindingTransformerWinding, nr::NetworkReductionData) =
    get(
        nr.impedance_correction_factors,
        (IS.get_uuid(get_transformer(w)), get_winding_number(w)),
        1.0,
    )

_impedance_corrected(br::PSY.ACTransmission, nr::NetworkReductionData) =
    !isone(_impedance_correction_factor(br, nr))
_impedance_corrected(bp::AbstractBranchesParallel, nr::NetworkReductionData) =
    any(br -> _impedance_corrected(br, nr), bp.branches)

# The DC model reads phase-shifting arcs and series-reduced chains from component
# reactances, which impedance correction does not touch; every other arc inherits the
# correction through Ybus. Warn about the mixed result until it is validated against PSS/E.
function _warn_impedance_correction_in_dc(nr::NetworkReductionData)
    any(!isone, values(nr.impedance_correction_factors)) || return
    corrected = Set(k[1] for (k, v) in nr.impedance_correction_factors if !isone(v))
    @warn "$(length(corrected)) transformer(s) carry an active impedance correction. " *
          "DC susceptances of phase-shifting transformers and series-reduced chains use " *
          "the uncorrected reactance; the DC model has not been validated against PSS/E " *
          "with impedance correction tables." maxlog = 1
    return
end

# DC susceptance of a symmetric (non-phase-shifting) arc from its Ybus off-diagonal, as
# `BA_Matrix` reads it: the reciprocal of the equivalent series reactance.
function _symmetric_arc_dc_susceptance(Y_ft::Complex)
    iszero(Y_ft) && return 0.0
    x_eq = imag(1 / Y_ft)
    iszero(x_eq) && return 0.0
    return 1 / x_eq
end

# The susceptance `BA_Matrix` assigns to an arc, from the same 2x2 it reads off Ybus. A
# phase-shifting arc (asymmetric off-diagonals) takes the phase-independent component value.
function _ba_arc_susceptance(entries::NTuple{4, <:Complex}, segment)
    Y_ft = -entries[2]
    Y_tf = -entries[3]
    Y_ft != Y_tf && return get_series_susceptance(segment)
    return _symmetric_arc_dc_susceptance(Y_ft)
end
