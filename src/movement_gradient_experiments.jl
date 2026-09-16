
# 1. Terminal N-R-C model:
#    Nutrient fluxes move directionally toward terminal node 5. Only N moves;
#    resources and consumers do not recycle back into the nutrient pool.
#
# 2. Well-mixed hub model:
#    Rosenzweig-MacArthur resource-consumer dynamics occur in each node, and
#    dispersal moves resources and consumers through a connected hub-satellite
#    network with directed satellite shortcuts.

using Statistics
using LinearAlgebra
using DelimitedFiles
using DifferentialEquations
using Random
using Plots

gr()

const TERMINAL_N_NODES = 5
const HUB_N_NODES = 5
const HUB = 1
const SATELLITES = [2, 3, 4, 5]
const SATELLITE_SHORTCUTS = [
    (2, 3),
    (3, 5),
    (4, 3),
    (5, 3),
]


const OUTPUT_DIR = normpath(joinpath(@__DIR__, "movement_gradient_output_final_management"))


const DISPLAY_PLOTS_IN_VSCODE = true

coefficient_of_variation(x; eps_value=1e-8) = begin
    m = mean(skipmissing(x))
    if !isfinite(m) || abs(m) <= eps_value
        return NaN
    end
    std(skipmissing(x)) / abs(m)
end

function save_table(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        writedlm(io, rows, ',')
    end
end

function node_outgoing_flux(x, edges_from, edges_to, node_rates)
    dx = zeros(length(x))
    for (source, target) in zip(edges_from, edges_to)
        flux = node_rates[source] * x[source]
        dx[source] -= flux
        dx[target] += flux
    end
    return dx
end

function adjacency_link_flux(x, adjacency, link_rates)
    dx = zeros(length(x))
    for source in axes(adjacency, 1), target in axes(adjacency, 2)
        if adjacency[source, target] != 0
            flux = link_rates[source, target] * x[source]
            dx[source] -= flux
            dx[target] += flux
        end
    end
    return dx
end

function post_transient_indices(times; transient=600.0)
    return findall(t -> t >= transient, times)
end

function percent_reduction_from_baseline(d_fast, d_slow)
    return 100.0 * (d_fast - d_slow) / d_fast
end

# ---------------------------------------------------------------------------
# Terminal N-R-C model

const TERMINAL_EDGES_FROM = [1, 2, 3, 4]
const TERMINAL_EDGES_TO = [2, 4, 4, 5]

struct NRCParams
    I::Float64
    loss::Float64
    uptake_rate::Float64
    half_saturation::Float64
    resource_mortality::Float64
    attack_rate::Float64
    handling_time::Float64
    conversion_efficiency::Float64
    consumer_mortality::Float64
    nutrient_flux_rates::Vector{Float64}
end

function terminal_nrc_rhs!(du, u, p::NRCParams, t)
    N = @view u[1:TERMINAL_N_NODES]
    R = @view u[(TERMINAL_N_NODES + 1):(2 * TERMINAL_N_NODES)]
    C = @view u[(2 * TERMINAL_N_NODES + 1):(3 * TERMINAL_N_NODES)]

    dN = @view du[1:TERMINAL_N_NODES]
    dR = @view du[(TERMINAL_N_NODES + 1):(2 * TERMINAL_N_NODES)]
    dC = @view du[(2 * TERMINAL_N_NODES + 1):(3 * TERMINAL_N_NODES)]

    for i in 1:TERMINAL_N_NODES
        Ni = max(N[i], eps(Float64))
        Ri = max(R[i], eps(Float64))
        Ci = max(C[i], eps(Float64))

        uptake = p.uptake_rate * (Ni / (p.half_saturation + Ni)) * Ri
        consumption = (p.attack_rate * Ri * Ci) /
            (1 + p.attack_rate * p.handling_time * Ri)

        dN[i] = p.I - p.loss * Ni - uptake
        dR[i] = uptake - p.resource_mortality * Ri - consumption
        dC[i] = p.conversion_efficiency * consumption - p.consumer_mortality * Ci
    end

    dN .+= node_outgoing_flux(N, TERMINAL_EDGES_FROM, TERMINAL_EDGES_TO, p.nutrient_flux_rates)
    return nothing
end

function make_nrc_params(node_rates)
    return NRCParams(
        0.4608592032838725,    # basal nutrient input
        0.015662473666718177,  # nutrient loss
        0.5349755765573964,    # uptake rate
        2.865132388775022,     # half saturation
        0.3203814387378674,    # resource mortality
        0.405356841253649,     # attack rate
        0.03904361181270938,   # handling time
        0.65,                  # conversion efficiency
        0.0757988537518309,    # consumer mortality
        copy(node_rates),
    )
end

function solve_terminal_nrc(node_rates; tspan=(0.0, 1200.0), saveat=0.05)
    u0 = vcat(
        fill(1.0, TERMINAL_N_NODES),
        fill(1.0, TERMINAL_N_NODES),
        fill(0.2, TERMINAL_N_NODES),
    )
    problem = ODEProblem(terminal_nrc_rhs!, u0, tspan, make_nrc_params(node_rates))
    solve(problem, Tsit5(); reltol=1e-8, abstol=1e-10, saveat=saveat, maxiters=10^7)
end

struct NRCPulseParams
    base_params::NRCParams
    base_flux_rate::Float64
    pulse_flux_addition::Float64
    pulse_period::Float64
    pulse_duration::Float64
    pulse_start::Float64
end

function annual_pulse_on(t, p)
    phase = mod(t - p.pulse_start, p.pulse_period)
    return phase < p.pulse_duration
end

function terminal_nrc_pulsed_rhs!(du, u, p::NRCPulseParams, t)
    flux_rate = p.base_flux_rate
    if annual_pulse_on(t, p)
        flux_rate += p.pulse_flux_addition
    end

    local_params = NRCParams(
        p.base_params.I,
        p.base_params.loss,
        p.base_params.uptake_rate,
        p.base_params.half_saturation,
        p.base_params.resource_mortality,
        p.base_params.attack_rate,
        p.base_params.handling_time,
        p.base_params.conversion_efficiency,
        p.base_params.consumer_mortality,
        fill(flux_rate, TERMINAL_N_NODES),
    )
    terminal_nrc_rhs!(du, u, local_params, t)
    return nothing
end

function annual_pulse_tstops(tspan, pulse_period, pulse_duration, pulse_start)
    stops = Float64[]
    pulse_begin = pulse_start
    while pulse_begin <= tspan[2]
        if pulse_begin >= tspan[1]
            push!(stops, pulse_begin)
        end
        pulse_end = pulse_begin + pulse_duration
        if pulse_end >= tspan[1] && pulse_end <= tspan[2]
            push!(stops, pulse_end)
        end
        pulse_begin += pulse_period
    end
    return sort(unique(stops))
end

function solve_terminal_nrc_annual_pulse(
    base_flux_rate;
    pulse_flux_addition=1.5,
    pulse_period=1.0,
    pulse_duration=0.1,
    pulse_start=0.0,
    tspan=(0.0, 1200.0),
    saveat=0.05,
)
    u0 = vcat(
        fill(1.0, TERMINAL_N_NODES),
        fill(1.0, TERMINAL_N_NODES),
        fill(0.2, TERMINAL_N_NODES),
    )
    params = NRCPulseParams(
        make_nrc_params(fill(base_flux_rate, TERMINAL_N_NODES)),
        base_flux_rate,
        pulse_flux_addition,
        pulse_period,
        pulse_duration,
        pulse_start,
    )
    problem = ODEProblem(terminal_nrc_pulsed_rhs!, u0, tspan, params)
    solve(
        problem,
        Tsit5();
        reltol=1e-8,
        abstol=1e-10,
        saveat=saveat,
        tstops=annual_pulse_tstops(tspan, pulse_period, pulse_duration, pulse_start),
        maxiters=10^7,
    )
end

function scan_terminal_nrc_speedup(; d_values=collect(0.0:0.05:1.25), transient=600.0)
    rows = Vector{NamedTuple}()
    for d in d_values
        sol = solve_terminal_nrc(fill(d, TERMINAL_N_NODES))
        indices = post_transient_indices(sol.t; transient)
        c5 = [sol.u[k][2 * TERMINAL_N_NODES + 5] for k in indices]
        push!(rows, (
            model="terminal_nrc",
            scenario="all_links",
            label="all nutrient links",
            movement_rate=d,
            metric=coefficient_of_variation(c5),
        ))
    end
    return rows
end

function scan_terminal_nrc_annual_pulses(;
    constant_flux_values=collect(0.0:0.05:1.25),
    pulse_multipliers=(4.0, 8.0),
    pulse_period=1.0,
    pulse_duration=0.1,
    pulse_start=0.0,
    transient=600.0,
)
    pulse_duty = pulse_duration / pulse_period
    treatments = [
        ("constant", "constant nutrient flux", 1.0, 0.0),
        ("pulse_4x", "4x pulse", pulse_multipliers[1], pulse_duration),
        ("pulse_8x", "8x pulse", pulse_multipliers[2], pulse_duration),
    ]

    rows = Vector{NamedTuple}()
    for constant_flux in constant_flux_values, (scenario, label, pulse_multiplier, active_duration) in treatments
        off_flux = constant_flux
        pulse_flux = pulse_multiplier * off_flux
        annual_mean_flux = if scenario == "constant"
            constant_flux
        else
            ((1.0 - pulse_duty) * off_flux) + (pulse_duty * pulse_flux)
        end

        sol = if scenario == "constant"
            solve_terminal_nrc(fill(constant_flux, TERMINAL_N_NODES))
        else
            solve_terminal_nrc_annual_pulse(
                off_flux;
                pulse_flux_addition=pulse_flux - off_flux,
                pulse_period=pulse_period,
                pulse_duration=pulse_duration,
                pulse_start=pulse_start,
            )
        end

        indices = post_transient_indices(sol.t; transient)
        n5 = [sol.u[k][5] for k in indices]
        r5 = [sol.u[k][TERMINAL_N_NODES + 5] for k in indices]
        c5 = [sol.u[k][2 * TERMINAL_N_NODES + 5] for k in indices]
        push!(rows, (
            model="terminal_nrc",
            scenario=scenario,
            label=label,
            constant_flux_rate=constant_flux,
            annual_mean_flux_rate=annual_mean_flux,
            off_pulse_flux_rate=off_flux,
            baseline_nutrient_input=make_nrc_params(fill(0.0, TERMINAL_N_NODES)).I,
            pulse_flux_multiplier=pulse_multiplier,
            pulse_flux_rate=pulse_flux,
            pulse_period=pulse_period,
            pulse_duration=active_duration,
            mean_N5=mean(n5),
            mean_R5=mean(r5),
            mean_C5=mean(c5),
            cv_C5=coefficient_of_variation(c5),
        ))
    end
    return rows
end

function scan_terminal_nrc_slowdown(; d_fast=1.25, d_values=reverse(collect(0.0:0.05:1.25)), transient=600.0)
    scenarios = [
        ("slow_link_4_to_5", "slow N link 4 -> 5"),
        ("mean_other_links", "mean across other N links"),
        ("all_links", "slow all N links"),
    ]

    rows = Vector{NamedTuple}()
    for (scenario, label) in scenarios, d_slow in d_values
        node_rates = fill(d_fast, TERMINAL_N_NODES)
        if scenario == "slow_link_4_to_5"
            node_rates[4] = d_slow
        elseif scenario == "all_links"
            node_rates .= d_slow
        end

        metric = if scenario == "mean_other_links"
            other_link_cvs = Float64[]
            for source in (1, 2, 3)
                single_link_rates = fill(d_fast, TERMINAL_N_NODES)
                single_link_rates[source] = d_slow
                sol = solve_terminal_nrc(single_link_rates)
                indices = post_transient_indices(sol.t; transient)
                c5 = [sol.u[k][2 * TERMINAL_N_NODES + 5] for k in indices]
                push!(other_link_cvs, coefficient_of_variation(c5))
            end
            mean(other_link_cvs)
        else
            sol = solve_terminal_nrc(node_rates)
            indices = post_transient_indices(sol.t; transient)
            c5 = [sol.u[k][2 * TERMINAL_N_NODES + 5] for k in indices]
            coefficient_of_variation(c5)
        end

        push!(rows, (
            model="terminal_nrc",
            scenario=scenario,
            label=label,
            d_slow=d_slow,
            percent_reduction=percent_reduction_from_baseline(d_fast, d_slow),
            metric=metric,
        ))
    end
    return rows
end

# ---------------------------------------------------------------------------
# Well-mixed hub-satellite Rosenzweig-MacArthur model

struct HubRMParams
    growth_rates::Vector{Float64}
    carrying_capacities::Vector{Float64}
    attack_rate::Float64
    handling_time::Float64
    conversion_efficiency::Float64
    mortality_rate::Float64
    resource_link_rates::Matrix{Float64}
    consumer_link_rates::Matrix{Float64}
    adjacency::Matrix{Int}
end

function well_mixed_hub_adjacency()
    adjacency = zeros(Int, HUB_N_NODES, HUB_N_NODES)
    for node in SATELLITES
        adjacency[HUB, node] = 1
        adjacency[node, HUB] = 1
    end
    for (source, target) in SATELLITE_SHORTCUTS
        adjacency[source, target] = 1
    end

    return adjacency, SATELLITE_SHORTCUTS
end

original_well_mixed_adjacency() = well_mixed_hub_adjacency()

function hub_rm_rhs!(du, u, p::HubRMParams, t)
    R = @view u[1:HUB_N_NODES]
    C = @view u[(HUB_N_NODES + 1):(2 * HUB_N_NODES)]
    dR = @view du[1:HUB_N_NODES]
    dC = @view du[(HUB_N_NODES + 1):(2 * HUB_N_NODES)]

    for i in 1:HUB_N_NODES
        Ri = max(R[i], eps(Float64))
        Ci = max(C[i], eps(Float64))
        consumption = (p.attack_rate * Ri * Ci) /
            (1 + p.attack_rate * p.handling_time * Ri)

        dR[i] = p.growth_rates[i] * Ri * (1 - Ri / p.carrying_capacities[i]) - consumption
        dC[i] = p.conversion_efficiency * consumption - p.mortality_rate * Ci
    end

    dR .+= adjacency_link_flux(R, p.adjacency, p.resource_link_rates)
    dC .+= adjacency_link_flux(C, p.adjacency, p.consumer_link_rates)
    return nothing
end

function link_rates_from_adjacency(adjacency, rate)
    return Float64.(adjacency) .* rate
end

function make_hub_rm_params(link_rates; adjacency)
    return HubRMParams(
        [
            3.1115294799548088,
            1.8500986097028593,
            1.5137170443023393,
            2.691052523204159,
            2.7751479145542888,
        ],
        [
            4.492631701792497,
            2.463701255821692,
            5.217249718210642,
            3.043395668956208,
            3.3332428755234655,
        ],
        0.7064665999286242, # attack rate
        1.963325466463706,  # handling time
        0.8,
        0.14878795596806443,
        copy(link_rates),
        copy(link_rates),
        adjacency,
    )
end

function solve_hub_rm(link_rates; adjacency, tspan=(0.0, 1600.0), saveat=0.25)
    params = make_hub_rm_params(link_rates; adjacency)
    u0 = vcat(0.8 .* params.carrying_capacities, 0.18 .* params.carrying_capacities)
    problem = ODEProblem(hub_rm_rhs!, u0, tspan, params)
    solve(problem, Tsit5(); reltol=1e-8, abstol=1e-10, saveat=saveat, maxiters=10^7)
end

struct HubRMPulseParams
    base_params::HubRMParams
    base_dispersal_rate::Float64
    pulse_dispersal_addition::Float64
    pulse_period::Float64
    pulse_duration::Float64
    pulse_start::Float64
end

function hub_rm_pulsed_rhs!(du, u, p::HubRMPulseParams, t)
    dispersal_rate = p.base_dispersal_rate
    if annual_pulse_on(t, p)
        dispersal_rate += p.pulse_dispersal_addition
    end

    link_rates = link_rates_from_adjacency(p.base_params.adjacency, dispersal_rate)
    local_params = HubRMParams(
        p.base_params.growth_rates,
        p.base_params.carrying_capacities,
        p.base_params.attack_rate,
        p.base_params.handling_time,
        p.base_params.conversion_efficiency,
        p.base_params.mortality_rate,
        link_rates,
        link_rates,
        p.base_params.adjacency,
    )
    hub_rm_rhs!(du, u, local_params, t)
    return nothing
end

function solve_hub_rm_annual_pulse(
    base_dispersal_rate;
    pulse_dispersal_addition=1.5,
    pulse_period=1.0,
    pulse_duration=0.1,
    pulse_start=0.0,
    adjacency,
    tspan=(0.0, 1600.0),
    saveat=0.25,
)
    base_params = make_hub_rm_params(
        link_rates_from_adjacency(adjacency, base_dispersal_rate);
        adjacency,
    )
    u0 = vcat(0.8 .* base_params.carrying_capacities, 0.18 .* base_params.carrying_capacities)
    params = HubRMPulseParams(
        base_params,
        base_dispersal_rate,
        pulse_dispersal_addition,
        pulse_period,
        pulse_duration,
        pulse_start,
    )
    problem = ODEProblem(hub_rm_pulsed_rhs!, u0, tspan, params)
    solve(
        problem,
        Tsit5();
        reltol=1e-8,
        abstol=1e-10,
        saveat=saveat,
        tstops=annual_pulse_tstops(tspan, pulse_period, pulse_duration, pulse_start),
        maxiters=10^7,
    )
end

function average_consumer_temporal_cv(sol; transient=800.0)
    indices = post_transient_indices(sol.t; transient)
    cvs = [
        coefficient_of_variation([sol.u[k][HUB_N_NODES + node] for k in indices])
        for node in 1:HUB_N_NODES
    ]
    return mean(skipmissing(cvs))
end

function mean_pairwise_consumer_synchrony(sol; transient=800.0)
    indices = post_transient_indices(sol.t; transient)
    pairwise_correlations = Float64[]

    for i in 1:(HUB_N_NODES - 1)
        for j in (i + 1):HUB_N_NODES
            x = [sol.u[k][HUB_N_NODES + i] for k in indices]
            y = [sol.u[k][HUB_N_NODES + j] for k in indices]

            if std(x) < eps(Float64) || std(y) < eps(Float64)
                push!(pairwise_correlations, NaN)
            else
                push!(pairwise_correlations, cor(x, y))
            end
        end
    end

    valid_correlations = filter(isfinite, pairwise_correlations)
    return isempty(valid_correlations) ? NaN : mean(valid_correlations)
end

function scan_hub_rm_speedup(; d_values=collect(0.0:0.05:1.25), transient=800.0)
    adjacency, satellite_shortcuts = well_mixed_hub_adjacency()
    rows = Vector{NamedTuple}()
    for d in d_values
        sol = solve_hub_rm(link_rates_from_adjacency(adjacency, d); adjacency)
        push!(rows, (
            model="hub_rm",
            scenario="all_links",
            label="all resource and consumer links",
            movement_rate=d,
            metric=average_consumer_temporal_cv(sol; transient),
            synchrony=mean_pairwise_consumer_synchrony(sol; transient),
        ))
    end
    return rows, adjacency, satellite_shortcuts
end

function scan_hub_rm_annual_pulses(;
    constant_dispersal_values=collect(0.0:0.05:1.25),
    pulse_multipliers=(4.0, 8.0),
    pulse_period=1.0,
    pulse_duration=0.1,
    pulse_start=0.0,
    transient=800.0,
)
    adjacency, satellite_shortcuts = well_mixed_hub_adjacency()
    pulse_duty = pulse_duration / pulse_period
    treatments = [
        ("constant", "constant dispersal", 1.0, 0.0),
        ("pulse_4x", "4x pulse", pulse_multipliers[1], pulse_duration),
        ("pulse_8x", "8x pulse", pulse_multipliers[2], pulse_duration),
    ]

    rows = Vector{NamedTuple}()
    for constant_dispersal in constant_dispersal_values, (scenario, label, pulse_multiplier, active_duration) in treatments
        off_dispersal = constant_dispersal
        pulse_dispersal = pulse_multiplier * off_dispersal
        annual_mean_dispersal = if scenario == "constant"
            constant_dispersal
        else
            ((1.0 - pulse_duty) * off_dispersal) + (pulse_duty * pulse_dispersal)
        end

        sol = if scenario == "constant"
            solve_hub_rm(link_rates_from_adjacency(adjacency, constant_dispersal); adjacency)
        else
            solve_hub_rm_annual_pulse(
                off_dispersal;
                pulse_dispersal_addition=pulse_dispersal - off_dispersal,
                pulse_period=pulse_period,
                pulse_duration=pulse_duration,
                pulse_start=pulse_start,
                adjacency,
            )
        end

        push!(rows, (
            model="hub_rm",
            scenario=scenario,
            label=label,
            constant_dispersal_rate=constant_dispersal,
            annual_mean_dispersal_rate=annual_mean_dispersal,
            off_pulse_dispersal_rate=off_dispersal,
            pulse_dispersal_multiplier=pulse_multiplier,
            pulse_dispersal_rate=pulse_dispersal,
            pulse_period=pulse_period,
            pulse_duration=active_duration,
            metric=average_consumer_temporal_cv(sol; transient),
            synchrony=mean_pairwise_consumer_synchrony(sol; transient),
        ))
    end
    return rows, adjacency, satellite_shortcuts
end

"""
Set the transfer rate for every directed currency flow entering or leaving
`node`. The model stores rates on directed links, so changing both incoming
and outgoing incident links is the explicit implementation of slowing
currency transfer at a node.
"""
function slow_currency_transfer_at_node(adjacency, d_fast, d_slow, node)
    transfer_rates = link_rates_from_adjacency(adjacency, d_fast)
    for source in axes(adjacency, 1), target in axes(adjacency, 2)
        if adjacency[source, target] != 0 && (source == node || target == node)
            transfer_rates[source, target] = d_slow
        end
    end
    return transfer_rates
end

"""Reduce transfer only on directed links leaving `node`."""
function slow_outgoing_currency_transfer_from_node(adjacency, d_fast, d_slow, node)
    transfer_rates = link_rates_from_adjacency(adjacency, d_fast)
    for target in axes(adjacency, 2)
        if adjacency[node, target] != 0
            transfer_rates[node, target] = d_slow
        end
    end
    return transfer_rates
end

function scan_hub_rm_slowdown(;
    d_fast=1.25,
    d_values=reverse(collect(0.0:0.05:1.25)),
    transient=800.0,
)
    adjacency, satellite_shortcuts = well_mixed_hub_adjacency()
    scenarios = [
        ("all_nodes", "slow currency transfer in/out of all nodes", "all"),
        ("hub_node", "slow currency transfer in/out of hub node 1", string(HUB)),
        (
            "mean_random_nodes",
            "mean in/out slowdown across satellite nodes",
            join(SATELLITES, ";"),
        ),
    ]

    rows = Vector{NamedTuple}()
    for (scenario, label, target_nodes) in scenarios, d_slow in d_values
        metric, synchrony = if scenario == "mean_random_nodes"
            node_metrics = Float64[]
            node_synchronies = Float64[]
            for node in SATELLITES
                transfer_rates = slow_currency_transfer_at_node(
                    adjacency,
                    d_fast,
                    d_slow,
                    node,
                )
                sol = solve_hub_rm(transfer_rates; adjacency)
                push!(node_metrics, average_consumer_temporal_cv(sol; transient))
                push!(node_synchronies, mean_pairwise_consumer_synchrony(sol; transient))
            end
            (mean(node_metrics), mean(node_synchronies))
        else
            transfer_rates = if scenario == "all_nodes"
                link_rates_from_adjacency(adjacency, d_slow)
            elseif scenario == "hub_node"
                slow_currency_transfer_at_node(adjacency, d_fast, d_slow, HUB)
            else
                error("Unknown node slowdown scenario: $scenario")
            end
            sol = solve_hub_rm(transfer_rates; adjacency)
            (
                average_consumer_temporal_cv(sol; transient),
                mean_pairwise_consumer_synchrony(sol; transient),
            )
        end

        push!(rows, (
            model="hub_rm",
            scenario=scenario,
            label=label,
            target_nodes=target_nodes,
            d_slow=d_slow,
            percent_reduction=percent_reduction_from_baseline(d_fast, d_slow),
            metric=metric,
            synchrony=synchrony,
        ))
    end
    return rows, adjacency, satellite_shortcuts
end

function validate_stabilize_then_destabilize(rows)
    for scenario in unique(row.scenario for row in rows)
        scenario_rows = sort(
            [row for row in rows if row.scenario == scenario && isfinite(row.metric)];
            by=row -> row.percent_reduction,
        )
        metrics = [row.metric for row in scenario_rows]
        minimum_index = argmin(metrics)
        baseline_metric = first(metrics)
        minimum_metric = metrics[minimum_index]
        final_metric = last(metrics)

        has_interior_minimum = 1 < minimum_index < length(metrics)
        initially_stabilizes = minimum_metric < baseline_metric
        subsequently_destabilizes = final_metric > minimum_metric
        if !(has_interior_minimum && initially_stabilizes && subsequently_destabilizes)
            error("Scenario $scenario does not stabilize and then destabilize")
        end

        println(
            "Validated ",
            scenario,
            ": CV ",
            round(baseline_metric; digits=4),
            " -> ",
            round(minimum_metric; digits=4),
            " -> ",
            round(final_metric; digits=4),
        )
    end
    return nothing
end

function slow_outgoing_links(adjacency, d_fast, d_slow, node)
    link_rates = link_rates_from_adjacency(adjacency, d_fast)
    for target in axes(adjacency, 2)
        if adjacency[node, target] != 0
            link_rates[node, target] = d_slow
        end
    end
    return link_rates
end

"""Reduce nutrient transfer only on directed links leaving `node`."""
function slow_terminal_outgoing_currency_transfer(d_fast, d_slow, node)
    node_rates = fill(d_fast, TERMINAL_N_NODES)
    node_rates[node] = d_slow
    return node_rates
end

function scan_terminal_nrc_reference_slowdown(;
    d_fast=1.25,
    d_values=reverse(collect(0.0:0.05:1.25)),
    transient=600.0,
)
    hub_node = 4
    non_hub_nonterminal_nodes = [1, 2, 3]
    scenarios = [
        ("all_nodes", "slow outgoing nutrient transfer from all nodes", "all"),
        ("hub_node", "slow outgoing nutrient transfer from hub node 4", string(hub_node)),
        (
            "mean_random_nodes",
            "mean outgoing slowdown across non-hub, non-terminal nodes",
            join(non_hub_nonterminal_nodes, ";"),
        ),
    ]

    rows = Vector{NamedTuple}()
    for (scenario, label, target_nodes) in scenarios, d_slow in d_values
        metric = if scenario == "mean_random_nodes"
            node_cvs = Float64[]
            for node in non_hub_nonterminal_nodes
                node_rates = slow_terminal_outgoing_currency_transfer(d_fast, d_slow, node)
                sol = solve_terminal_nrc(node_rates)
                indices = post_transient_indices(sol.t; transient)
                c5 = [sol.u[k][2 * TERMINAL_N_NODES + 5] for k in indices]
                push!(node_cvs, coefficient_of_variation(c5))
            end
            mean(node_cvs)
        else
            node_rates = if scenario == "all_nodes"
                fill(d_slow, TERMINAL_N_NODES)
            elseif scenario == "hub_node"
                slow_terminal_outgoing_currency_transfer(d_fast, d_slow, hub_node)
            else
                error("Unknown terminal node slowdown scenario: $scenario")
            end
            sol = solve_terminal_nrc(node_rates)
            indices = post_transient_indices(sol.t; transient)
            c5 = [sol.u[k][2 * TERMINAL_N_NODES + 5] for k in indices]
            coefficient_of_variation(c5)
        end
        push!(rows, (
            model="terminal_nrc",
            scenario=scenario,
            label=label,
            target_nodes=target_nodes,
            d_slow=d_slow,
            percent_reduction=percent_reduction_from_baseline(d_fast, d_slow),
            metric=metric,
        ))
    end
    return rows
end

function scan_hub_rm_reference_slowdown(; d_fast=1.25, d_values=reverse(collect(0.0:0.05:1.25)), transient=800.0)
    adjacency, satellite_shortcuts = well_mixed_hub_adjacency()
    scenarios = [
        ("slow_hub_1", "slow all links to/from hub"),
        ("slow_node_3", "slow all links to/from node 3"),
        ("slow_from_node_3", "slow all links from node 3"),
    ]

    rows = Vector{NamedTuple}()
    for (scenario, label) in scenarios, d_slow in d_values
        link_rates = if scenario == "slow_hub_1"
            slow_currency_transfer_at_node(adjacency, d_fast, d_slow, HUB)
        elseif scenario == "slow_node_3"
            slow_currency_transfer_at_node(adjacency, d_fast, d_slow, 3)
        elseif scenario == "slow_from_node_3"
            slow_outgoing_links(adjacency, d_fast, d_slow, 3)
        else
            link_rates_from_adjacency(adjacency, d_fast)
        end

        sol = solve_hub_rm(link_rates; adjacency)
        push!(rows, (
            model="hub_rm",
            scenario=scenario,
            label=label,
            d_slow=d_slow,
            percent_reduction=percent_reduction_from_baseline(d_fast, d_slow),
            metric=average_consumer_temporal_cv(sol; transient),
            synchrony=mean_pairwise_consumer_synchrony(sol; transient),
        ))
    end
    return rows, adjacency, satellite_shortcuts
end

# ---------------------------------------------------------------------------
# Output

const COLORS = Dict(
    "all_nodes" => "#5F5F5F",
    "all_links" => "#5F5F5F",
    "slow_link_4_to_5" => "#3568B8",
    "slow_link_3_to_4" => "#D9271C",
    "mean_random_links" => "#D9271C",
    "slow_hub_1" => "#3568B8",
    "slow_node_3" => "#D9271C",
    "hub_node" => "#3568B8",
    "mean_random_nodes" => "#D9271C",
    "slow_from_node_3" => "#D9271C",
    "slow_all_hub_links" => "#5F5F5F",
    "constant" => "#5F5F5F",
    "pulse_4x" => "#D99A32",
    "pulse_8x" => "#C83F2E",
)

const FONT_FAMILY = "Helvetica"

const LINE_STYLES = Dict(
    "all_links" => :solid,
    "slow_link_4_to_5" => :solid,
    "slow_link_3_to_4" => :solid,
    "mean_random_links" => :solid,
    "slow_hub_1" => :solid,
    "slow_node_3" => :solid,
    "all_nodes" => :solid,
    "hub_node" => :solid,
    "mean_random_nodes" => :solid,
    "slow_from_node_3" => :solid,
    "slow_all_hub_links" => :solid,
    "constant" => :solid,
    "pulse_4x" => :solid,
    "pulse_8x" => :solid,
)

function rows_to_matrix(rows)
    matrix = Matrix{Any}(undef, length(rows), 5)
    for (i, row) in enumerate(rows)
        x = hasproperty(row, :d_slow) ? row.d_slow : row.movement_rate
        matrix[i, :] .= (row.model, row.scenario, row.label, x, row.metric)
    end
    return matrix
end

function mitigation_rows_to_matrix(rows)
    matrix = Matrix{Any}(undef, length(rows), 6)
    for (i, row) in enumerate(rows)
        matrix[i, :] .= (
            row.model,
            row.scenario,
            row.label,
            row.d_slow,
            row.percent_reduction,
            row.metric,
        )
    end
    return matrix
end

function synchrony_rows_to_matrix(rows)
    matrix = Matrix{Any}(undef, length(rows), 6)
    for (i, row) in enumerate(rows)
        matrix[i, :] .= (
            row.model,
            row.scenario,
            row.label,
            row.d_slow,
            row.percent_reduction,
            row.synchrony,
        )
    end
    return matrix
end

function node_slowdown_rows_to_matrix(rows; response_field=:metric)
    matrix = Matrix{Any}(undef, length(rows), 7)
    for (i, row) in enumerate(rows)
        response = response_field == :metric ? row.metric : row.synchrony
        matrix[i, :] .= (
            row.model,
            row.scenario,
            row.label,
            row.target_nodes,
            row.d_slow,
            row.percent_reduction,
            response,
        )
    end
    return matrix
end

function speedup_synchrony_rows_to_matrix(rows)
    matrix = Matrix{Any}(undef, length(rows), 5)
    for (i, row) in enumerate(rows)
        matrix[i, :] .= (
            row.model,
            row.scenario,
            row.label,
            row.movement_rate,
            row.synchrony,
        )
    end
    return matrix
end

function pulse_rows_to_matrix(rows)
    matrix = Matrix{Any}(undef, length(rows), 14)
    for (i, row) in enumerate(rows)
        matrix[i, :] .= (
            row.model,
            row.scenario,
            row.label,
            row.constant_flux_rate,
            row.annual_mean_flux_rate,
            row.off_pulse_flux_rate,
            row.baseline_nutrient_input,
            row.pulse_flux_multiplier,
            row.pulse_flux_rate,
            row.pulse_period,
            row.pulse_duration,
            row.mean_N5,
            row.mean_R5,
            row.mean_C5,
        )
    end
    return matrix
end

function pulse_cv_rows_to_matrix(rows)
    matrix = Matrix{Any}(undef, length(rows), 15)
    for (i, row) in enumerate(rows)
        matrix[i, :] .= (
            row.model,
            row.scenario,
            row.label,
            row.constant_flux_rate,
            row.annual_mean_flux_rate,
            row.off_pulse_flux_rate,
            row.baseline_nutrient_input,
            row.pulse_flux_multiplier,
            row.pulse_flux_rate,
            row.pulse_period,
            row.pulse_duration,
            row.mean_N5,
            row.mean_R5,
            row.mean_C5,
            row.cv_C5,
        )
    end
    return matrix
end

function hub_pulse_rows_to_matrix(rows)
    matrix = Matrix{Any}(undef, length(rows), 12)
    for (i, row) in enumerate(rows)
        matrix[i, :] .= (
            row.model,
            row.scenario,
            row.label,
            row.constant_dispersal_rate,
            row.annual_mean_dispersal_rate,
            row.off_pulse_dispersal_rate,
            row.pulse_dispersal_multiplier,
            row.pulse_dispersal_rate,
            row.pulse_period,
            row.pulse_duration,
            row.metric,
            row.synchrony,
        )
    end
    return matrix
end

function save_terminal_network(output_path)
    rows = hcat(TERMINAL_EDGES_FROM, TERMINAL_EDGES_TO)
    save_table(output_path, ["source", "target"], rows)
end

function save_well_mixed_network(adjacency, satellite_shortcuts, output_path)
    rows = Matrix{Any}(undef, count(x -> x != 0, adjacency), 3)
    shortcut_edges = Set(satellite_shortcuts)
    row = 1
    for source in axes(adjacency, 1), target in axes(adjacency, 2)
        if adjacency[source, target] != 0
            link_type = if source == HUB || target == HUB
                "hub_primary"
            elseif (source, target) in shortcut_edges
                "satellite_shortcut"
            else
                "satellite_link"
            end
            rows[row, :] .= (source, target, link_type)
            row += 1
        end
    end
    save_table(output_path, ["source", "target", "link_type"], rows)
end

function save_parameter_table(output_path)
    nrc = make_nrc_params(fill(0.0, TERMINAL_N_NODES))
    nrc_rows = [
        ("terminal_nrc", "basal_nutrient_input", nrc.I),
        ("terminal_nrc", "nutrient_loss", nrc.loss),
        ("terminal_nrc", "uptake_rate", nrc.uptake_rate),
        ("terminal_nrc", "half_saturation", nrc.half_saturation),
        ("terminal_nrc", "resource_mortality", nrc.resource_mortality),
        ("terminal_nrc", "attack_rate", nrc.attack_rate),
        ("terminal_nrc", "handling_time", nrc.handling_time),
        ("terminal_nrc", "conversion_efficiency", nrc.conversion_efficiency),
        ("terminal_nrc", "consumer_mortality", nrc.consumer_mortality),
    ]

    adjacency, _ = well_mixed_hub_adjacency()
    hub = make_hub_rm_params(link_rates_from_adjacency(adjacency, 0.0); adjacency)
    hub_rows = [
        ("hub_cr_metacommunity", "attack_rate", hub.attack_rate),
        ("hub_cr_metacommunity", "handling_time", hub.handling_time),
        ("hub_cr_metacommunity", "conversion_efficiency", hub.conversion_efficiency),
        ("hub_cr_metacommunity", "mortality_rate", hub.mortality_rate),
    ]
    for node in 1:HUB_N_NODES
        push!(hub_rows, ("hub_cr_metacommunity", "growth_rate_node_$node", hub.growth_rates[node]))
        push!(hub_rows, ("hub_cr_metacommunity", "carrying_capacity_node_$node", hub.carrying_capacities[node]))
    end

    rows = vcat(nrc_rows, hub_rows)
    matrix = Matrix{Any}(undef, length(rows), 3)
    for (i, row) in enumerate(rows)
        matrix[i, :] .= row
    end
    save_table(output_path, ["model", "parameter", "value"], matrix)
end

function terminal_node_coordinates()
    return Dict(
        1 => (-1.05, 0.34),
        2 => (-0.46, 0.34),
        3 => (-0.46, -0.34),
        4 => (0.34, 0.0),
        5 => (1.0, 0.0),
    )
end

function hub_node_coordinates()
    return Dict(
        1 => (0.0, 0.0),
        2 => (-0.86, 0.52),
        3 => (0.86, 0.52),
        4 => (0.86, -0.52),
        5 => (-0.86, -0.52),
    )
end

function draw_directed_edge!(p, coords, source, target; color="#444444", linewidth=2.5, linestyle=:solid, alpha=1.0)
    x1, y1 = coords[source]
    x2, y2 = coords[target]
    dx = x2 - x1
    dy = y2 - y1
    edge_length = sqrt(dx^2 + dy^2)
    ux = dx / edge_length
    uy = dy / edge_length
    node_radius = 0.13
    tip_x = x2 - node_radius * ux
    tip_y = y2 - node_radius * uy
    base_x = x1 + node_radius * ux
    base_y = y1 + node_radius * uy
    plot!(
        p,
        [base_x, tip_x],
        [base_y, tip_y];
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        alpha=alpha,
        label=false,
    )
    head_length = 0.12
    head_width = 0.07
    left_x = tip_x - head_length * ux - head_width * uy
    left_y = tip_y - head_length * uy + head_width * ux
    right_x = tip_x - head_length * ux + head_width * uy
    right_y = tip_y - head_length * uy - head_width * ux
    plot!(p, [left_x, tip_x, right_x], [left_y, tip_y, right_y];
        color=color, linewidth=linewidth, linestyle=linestyle, alpha=alpha, label=false)
end

function draw_undirected_edge!(p, coords, source, target; color="#444444", linewidth=2.5, linestyle=:solid, alpha=1.0)
    x1, y1 = coords[source]
    x2, y2 = coords[target]
    plot!(
        p,
        [x1, x2],
        [y1, y2];
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        alpha=alpha,
        label=false,
    )
end

function plot_network_configurations(output_stem)
    mkpath(dirname(output_stem))
    terminal_coords = terminal_node_coordinates()
    hub_coords = hub_node_coordinates()
    hub_adjacency, satellite_shortcuts = well_mixed_hub_adjacency()

    terminal_plot = plot(
        title="A. Directionally-skewed network",
        xlim=(-1.25, 1.25),
        ylim=(-1.05, 1.05),
        aspect_ratio=:equal,
        axis=false,
        grid=false,
        legend=false,
        size=(900, 420),
        titlefontsize=15,
        background_color=:white,
    )
    for (source, target) in zip(TERMINAL_EDGES_FROM, TERMINAL_EDGES_TO)
        draw_directed_edge!(terminal_plot, terminal_coords, source, target; color="#111111", linewidth=2.6)
    end
    for node in 1:TERMINAL_N_NODES
        x, y = terminal_coords[node]
        markersize = node == 5 ? 22 : 15
        scatter!(terminal_plot, [x], [y]; markersize=markersize, markercolor="#000000",
            markerstrokecolor="#000000", markerstrokewidth=1, label=false)
    end

    hub_plot = plot(
        title="B. Well-mixed network",
        xlim=(-1.25, 1.25),
        ylim=(-1.25, 1.25),
        aspect_ratio=:equal,
        axis=false,
        grid=false,
        legend=false,
        size=(900, 420),
        titlefontsize=15,
        background_color=:white,
    )
    for node in SATELLITES
        draw_directed_edge!(hub_plot, hub_coords, HUB, node; color="#9E9E9E", linewidth=1.7, alpha=0.85)
        draw_directed_edge!(hub_plot, hub_coords, node, HUB; color="#9E9E9E", linewidth=1.7, alpha=0.85)
    end
    for (source, target) in satellite_shortcuts
        draw_directed_edge!(
            hub_plot,
            hub_coords,
            source,
            target;
            color="#111111",
            linewidth=1.8,
        )
    end
    for node in 1:HUB_N_NODES
        x, y = hub_coords[node]
        markersize = node == HUB ? 25 : 15
        scatter!(hub_plot, [x], [y]; markersize=markersize, markercolor="#000000",
            markerstrokecolor="#000000", markerstrokewidth=1, label=false)
    end

    combined = plot(terminal_plot, hub_plot; layout=(2, 1), size=(900, 840))
    png_path = output_stem * ".png"
    pdf_path = output_stem * ".pdf"
    savefig(combined, png_path)
    savefig(combined, pdf_path)
    if DISPLAY_PLOTS_IN_VSCODE
        display(combined)
    end
    println("Saved plot: " * png_path)
    println("Saved plot: " * pdf_path)
    return combined
end

function summarize_anchor_values(rows, movement_rates)
    for d in movement_rates
        closest = rows[argmin(abs.(Float64[row.movement_rate for row in rows] .- d))]
        println(
            closest.model,
            " CV at movement rate ",
            closest.movement_rate,
            " = ",
            round(closest.metric; digits=4),
        )
    end
end

function row_x_value(row, xfield)
    if xfield == :movement_rate
        return row.movement_rate
    elseif xfield == :dominant_flux_rate
        return row.dominant_flux_rate
    elseif xfield == :annual_mean_flux_rate
        return row.annual_mean_flux_rate
    elseif xfield == :constant_flux_rate
        return row.constant_flux_rate
    elseif xfield == :constant_dispersal_rate
        return row.constant_dispersal_rate
    elseif xfield == :percent_reduction
        return row.percent_reduction
    elseif xfield == :d_slow
        return row.d_slow
    end
    error("Unsupported xfield: $xfield")
end

function row_y_value(row, yfield)
    if yfield == :metric
        return row.metric
    elseif yfield == :synchrony
        return row.synchrony
    elseif yfield == :mean_C5
        return row.mean_C5
    elseif yfield == :cv_C5
        return row.cv_C5
    end
    error("Unsupported yfield: " * string(yfield))
end

function plot_rows(rows, output_stem; title="", xlabel, ylabel, xfield=:movement_rate, yfield=:metric, legend=false)
    mkpath(dirname(output_stem))

    scenario_order = unique([row.scenario for row in rows])
    p = plot(
        xlabel=xlabel,
        ylabel=ylabel,
        title="",
        legend=false,
        grid=:y,
        gridcolor="#E7E2DA",
        gridalpha=0.36,
        framestyle=:axes,
        background_color=:white,
        foreground_color=:black,
        foreground_color_text=:black,
        foreground_color_axis=:black,
        foreground_color_border=:black,
        size=(820, 560),
        fontfamily=FONT_FAMILY,
        guidefontsize=18,
        tickfontsize=14,
        legendfontsize=14,
        guidefontcolor=:black,
        tickfontcolor=:black,
        legendfontcolor=:black,
        legend_background_color=:white,
        legend_foreground_color=:black,
        legend_frame=false,
        left_margin=7Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=2Plots.mm,
        right_margin=4Plots.mm,
    )

    for scenario in scenario_order
        subset = [row for row in rows if row.scenario == scenario]
        subset = sort(subset, by=row -> row_x_value(row, xfield))
        line_style = length(scenario_order) == 1 ? :solid : get(LINE_STYLES, scenario, :solid)
        plot!(
            p,
            [row_x_value(row, xfield) for row in subset],
            [row_y_value(row, yfield) for row in subset];
            label=false,
            linewidth=5.4,
            color=COLORS[scenario],
            linealpha=0.88,
            linestyle=line_style,
            marker=:none,
        )
    end

    pdf_path = output_stem * ".pdf"
    png_path = output_stem * ".png"
    savefig(p, pdf_path)
    savefig(p, png_path)

    if !isfile(png_path)
        error("PNG was not saved: " * png_path)
    end

    if DISPLAY_PLOTS_IN_VSCODE
        display(p)
    end

    println("Saved plot: " * png_path)
    println("Saved plot: " * pdf_path)
    return p
end

function plot_terminal_pulse_figure_s1(rows, output_stem)
    mkpath(dirname(output_stem))

    label_text = Dict(
        "constant" => "constant",
        "pulse_4x" => "4x pulse",
        "pulse_8x" => "8x pulse",
    )

    cv_label_y_positions = Dict(
        "constant" => 0.342,
        "pulse_4x" => 0.413,
        "pulse_8x" => 0.475,
    )
    cv_panel = plot(
        xlabel="Directional nutrient flux rate",
        ylabel="Node 5 consumer CV",
        title="A",
        xlim=(0.0, 1.55),
        legend=false,
        grid=:y,
        gridcolor="#E7E2DA",
        gridalpha=0.36,
        framestyle=:axes,
        background_color=:white,
        foreground_color=:black,
        size=(820, 560),
        fontfamily=FONT_FAMILY,
        guidefontsize=15,
        tickfontsize=12,
        titlefontsize=16,
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=3Plots.mm,
        right_margin=7Plots.mm,
    )

    for scenario in ("constant", "pulse_4x", "pulse_8x")
        subset = sort([row for row in rows if row.scenario == scenario], by=row -> row.constant_flux_rate)
        plot!(
            cv_panel,
            [row.constant_flux_rate for row in subset],
            [row.cv_C5 for row in subset];
            color=COLORS[scenario],
            linewidth=4.8,
            linealpha=0.88,
            label=false,
        )
        annotate!(
            cv_panel,
            1.29,
            cv_label_y_positions[scenario],
            text("  " * label_text[scenario], 10, :left),
        )
    end

    constant_by_flux = Dict(
        row.constant_flux_rate => row.cv_C5
        for row in rows if row.scenario == "constant"
    )
    delta_panel = plot(
        xlabel="Directional nutrient flux rate",
        ylabel="Change in node 5 consumer CV",
        title="B",
        xlim=(0.0, 1.55),
        legend=false,
        grid=:y,
        gridcolor="#E7E2DA",
        gridalpha=0.36,
        framestyle=:axes,
        background_color=:white,
        foreground_color=:black,
        size=(820, 560),
        fontfamily=FONT_FAMILY,
        guidefontsize=15,
        tickfontsize=12,
        titlefontsize=16,
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=3Plots.mm,
        right_margin=7Plots.mm,
    )
    hline!(delta_panel, [0.0]; color="#5F5F5F", linewidth=2.0, linealpha=0.55, label=false)

    for scenario in ("pulse_4x", "pulse_8x")
        subset = sort([row for row in rows if row.scenario == scenario], by=row -> row.constant_flux_rate)
        deltas = [row.cv_C5 - constant_by_flux[row.constant_flux_rate] for row in subset]
        plot!(
            delta_panel,
            [row.constant_flux_rate for row in subset],
            deltas;
            color=COLORS[scenario],
            linewidth=4.8,
            linealpha=0.88,
            label=false,
        )
        annotate!(
            delta_panel,
            subset[end].constant_flux_rate,
            deltas[end],
            text("  " * label_text[scenario], 10, :left),
        )
    end

    figure = plot(cv_panel, delta_panel; layout=(1, 2), size=(1500, 560))
    png_path = output_stem * ".png"
    pdf_path = output_stem * ".pdf"
    savefig(figure, png_path)
    savefig(figure, pdf_path)
    if DISPLAY_PLOTS_IN_VSCODE
        display(figure)
    end
    println("Saved plot: " * png_path)
    println("Saved plot: " * pdf_path)
    return figure
end

function plot_hub_pulse_figure(rows, output_stem)
    mkpath(dirname(output_stem))

    label_text = Dict(
        "constant" => "constant",
        "pulse_4x" => "4x pulse",
        "pulse_8x" => "8x pulse",
    )

    cv_panel = plot(
        xlabel="Constant dispersal rate",
        ylabel="Average consumer CV",
        title="A",
        xlim=(0.0, 1.55),
        legend=false,
        grid=:y,
        gridcolor="#E7E2DA",
        gridalpha=0.36,
        framestyle=:axes,
        background_color=:white,
        foreground_color=:black,
        size=(820, 560),
        fontfamily=FONT_FAMILY,
        guidefontsize=15,
        tickfontsize=12,
        titlefontsize=16,
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=3Plots.mm,
        right_margin=7Plots.mm,
    )

    for scenario in ("constant", "pulse_4x", "pulse_8x")
        subset = sort([row for row in rows if row.scenario == scenario], by=row -> row.constant_dispersal_rate)
        x = [row.constant_dispersal_rate for row in subset]
        y = [row.metric for row in subset]
        plot!(
            cv_panel,
            x,
            y;
            color=COLORS[scenario],
            linewidth=4.8,
            linealpha=0.88,
            label=false,
        )
        annotate!(
            cv_panel,
            x[end],
            y[end],
            text("  " * label_text[scenario], 10, :left),
        )
    end

    constant_by_dispersal = Dict(
        row.constant_dispersal_rate => row.metric
        for row in rows if row.scenario == "constant"
    )
    delta_panel = plot(
        xlabel="Constant dispersal rate",
        ylabel="Change in average consumer CV",
        title="B",
        xlim=(0.0, 1.55),
        legend=false,
        grid=:y,
        gridcolor="#E7E2DA",
        gridalpha=0.36,
        framestyle=:axes,
        background_color=:white,
        foreground_color=:black,
        size=(820, 560),
        fontfamily=FONT_FAMILY,
        guidefontsize=15,
        tickfontsize=12,
        titlefontsize=16,
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=3Plots.mm,
        right_margin=7Plots.mm,
    )
    hline!(delta_panel, [0.0]; color="#5F5F5F", linewidth=2.0, linealpha=0.55, label=false)

    for scenario in ("pulse_4x", "pulse_8x")
        subset = sort([row for row in rows if row.scenario == scenario], by=row -> row.constant_dispersal_rate)
        x = [row.constant_dispersal_rate for row in subset]
        deltas = [row.metric - constant_by_dispersal[row.constant_dispersal_rate] for row in subset]
        plot!(
            delta_panel,
            x,
            deltas;
            color=COLORS[scenario],
            linewidth=4.8,
            linealpha=0.88,
            label=false,
        )
        annotate!(
            delta_panel,
            x[end],
            deltas[end],
            text("  " * label_text[scenario], 10, :left),
        )
    end

    figure = plot(cv_panel, delta_panel; layout=(1, 2), size=(1500, 560))
    png_path = output_stem * ".png"
    pdf_path = output_stem * ".pdf"
    savefig(figure, png_path)
    savefig(figure, pdf_path)
    if DISPLAY_PLOTS_IN_VSCODE
        display(figure)
    end
    println("Saved plot: " * png_path)
    println("Saved plot: " * pdf_path)
    return figure
end

function main()
    mkpath(OUTPUT_DIR)
    println("Running final hybrid management implementation")
    println("Directional model: outgoing transfer only")
    println("Well-mixed model: incoming and outgoing transfer")
    println("Slowdown scenarios: all_nodes, hub_node, mean_random_nodes")
    println("Fresh outputs will be written to " * OUTPUT_DIR)

    terminal_speedup_rows = scan_terminal_nrc_speedup()
    terminal_slowdown_rows = scan_terminal_nrc_reference_slowdown()
    terminal_pulse_rows = scan_terminal_nrc_annual_pulses()
    hub_speedup_rows, hub_adjacency, hub_satellite_choices = scan_hub_rm_speedup()
    hub_slowdown_rows, _, _ = scan_hub_rm_slowdown()
    validate_stabilize_then_destabilize(hub_slowdown_rows)
    hub_pulse_rows, _, _ = scan_hub_rm_annual_pulses()

    save_table(
        joinpath(OUTPUT_DIR, "terminal_nrc_speedup_scan.csv"),
        ["model", "scenario", "label", "movement_rate", "cv_C5"],
        rows_to_matrix(terminal_speedup_rows),
    )

    save_table(
        joinpath(OUTPUT_DIR, "terminal_nrc_slowdown_scan.csv"),
        [
            "model",
            "scenario",
            "label",
            "target_nodes",
            "d_slow",
            "percent_reduction",
            "cv_C5",
        ],
        node_slowdown_rows_to_matrix(terminal_slowdown_rows),
    )

    save_table(
        joinpath(OUTPUT_DIR, "terminal_nrc_annual_connectivity_pulse_scan.csv"),
        [
            "model",
            "scenario",
            "label",
            "constant_flux_rate",
            "annual_mean_flux_rate",
            "off_pulse_flux_rate",
            "baseline_nutrient_input",
            "pulse_flux_multiplier",
            "pulse_flux_rate",
            "pulse_period",
            "pulse_duration",
            "mean_N5",
            "mean_R5",
            "mean_C5",
            "cv_C5",
        ],
        pulse_cv_rows_to_matrix(terminal_pulse_rows),
    )

    save_table(
        joinpath(OUTPUT_DIR, "hub_rm_speedup_scan.csv"),
        ["model", "scenario", "label", "movement_rate", "average_consumer_temporal_cv"],
        rows_to_matrix(hub_speedup_rows),
    )

    save_table(
        joinpath(OUTPUT_DIR, "hub_rm_speedup_synchrony_scan.csv"),
        ["model", "scenario", "label", "movement_rate", "average_pairwise_consumer_synchrony"],
        speedup_synchrony_rows_to_matrix(hub_speedup_rows),
    )

    save_table(
        joinpath(OUTPUT_DIR, "hub_rm_slowdown_scan.csv"),
        [
            "model",
            "scenario",
            "label",
            "target_nodes",
            "d_slow",
            "percent_reduction",
            "average_consumer_temporal_cv",
        ],
        node_slowdown_rows_to_matrix(hub_slowdown_rows),
    )

    save_table(
        joinpath(OUTPUT_DIR, "hub_rm_slowdown_synchrony_scan.csv"),
        [
            "model",
            "scenario",
            "label",
            "target_nodes",
            "d_slow",
            "percent_reduction",
            "average_pairwise_consumer_synchrony",
        ],
        node_slowdown_rows_to_matrix(hub_slowdown_rows; response_field=:synchrony),
    )

    save_table(
        joinpath(OUTPUT_DIR, "hub_rm_annual_connectivity_pulse_scan.csv"),
        [
            "model",
            "scenario",
            "label",
            "constant_dispersal_rate",
            "annual_mean_dispersal_rate",
            "off_pulse_dispersal_rate",
            "pulse_dispersal_multiplier",
            "pulse_dispersal_rate",
            "pulse_period",
            "pulse_duration",
            "average_consumer_temporal_cv",
            "average_pairwise_consumer_synchrony",
        ],
        hub_pulse_rows_to_matrix(hub_pulse_rows),
    )

    save_terminal_network(joinpath(OUTPUT_DIR, "terminal_nrc_network_edges.csv"))
    save_well_mixed_network(
        hub_adjacency,
        hub_satellite_choices,
        joinpath(OUTPUT_DIR, "hub_rm_network_edges.csv"),
    )
    save_parameter_table(joinpath(OUTPUT_DIR, "selected_parameters.csv"))
    plot_network_configurations(joinpath(OUTPUT_DIR, "network_configurations"))

    plot_rows(
        terminal_speedup_rows,
        joinpath(OUTPUT_DIR, "terminal_nrc_speedup_scan");
        title="Directional terminal N-R-C model: nutrient movement gradient",
        xlabel="Directional nutrient flux rate",
        ylabel="CV of C5",
        legend=false,
    )
    plot_rows(
        terminal_slowdown_rows,
        joinpath(OUTPUT_DIR, "terminal_nrc_slowdown_scan");
        title="Directional terminal N-R-C model: targeted nutrient slow-down",
        xlabel="Reduction in outgoing nutrient-transfer rate (%)",
        ylabel="Consumer CV",
        xfield=:percent_reduction,
        legend=:bottomleft,
    )
    plot_rows(
        terminal_pulse_rows,
        joinpath(OUTPUT_DIR, "terminal_nrc_annual_connectivity_pulse_mean_scan");
        xlabel="Directional nutrient flux rate",
        ylabel="Mean C5",
        xfield=:constant_flux_rate,
        yfield=:mean_C5,
        legend=false,
    )
    plot_rows(
        terminal_pulse_rows,
        joinpath(OUTPUT_DIR, "terminal_nrc_annual_connectivity_pulse_cv_scan");
        xlabel="Directional nutrient flux rate",
        ylabel="CV of C5",
        xfield=:constant_flux_rate,
        yfield=:cv_C5,
        legend=false,
    )
    plot_terminal_pulse_figure_s1(
        terminal_pulse_rows,
        joinpath(OUTPUT_DIR, "figure_s1_temporal_pulse_amplification"),
    )
    plot_rows(
        hub_pulse_rows,
        joinpath(OUTPUT_DIR, "hub_rm_annual_connectivity_pulse_cv_scan");
        xlabel="Constant dispersal rate",
        ylabel="Average consumer CV",
        xfield=:constant_dispersal_rate,
        legend=false,
    )
    plot_hub_pulse_figure(
        hub_pulse_rows,
        joinpath(OUTPUT_DIR, "hub_rm_annual_connectivity_pulse_amplification"),
    )
    plot_rows(
        hub_speedup_rows,
        joinpath(OUTPUT_DIR, "hub_rm_speedup_scan");
        title="Well-mixed hub CR metacommunity: dispersal gradient",
        xlabel="Resource and consumer dispersal rate across all links",
        ylabel="Average consumer temporal CV",
        legend=false,
    )
    plot_rows(
        hub_speedup_rows,
        joinpath(OUTPUT_DIR, "hub_rm_speedup_synchrony_scan");
        xlabel="Resource and consumer dispersal rate across all links",
        ylabel="Average pairwise consumer synchrony",
        yfield=:synchrony,
        legend=false,
    )
    plot_rows(
        hub_slowdown_rows,
        joinpath(OUTPUT_DIR, "hub_rm_slowdown_scan");
        title="Well-mixed hub CR metacommunity: node currency-transfer slowdown",
        xlabel="Reduction in node currency-transfer rate (%)",
        ylabel="Average consumer CV",
        xfield=:percent_reduction,
        legend=:topleft,
    )
    plot_rows(
        hub_slowdown_rows,
        joinpath(OUTPUT_DIR, "hub_rm_slowdown_synchrony_scan");
        xlabel="Reduction in node currency-transfer rate (%)",
        ylabel="Average pairwise consumer synchrony",
        xfield=:percent_reduction,
        yfield=:synchrony,
        legend=:outerright,
    )

    println()
    println("Anchor checks for the desired checkmark CV response:")
    summarize_anchor_values(terminal_speedup_rows, (0.0, 0.45, 1.25))
    summarize_anchor_values(hub_speedup_rows, (0.0, 0.45, 1.25))
    println()
    println("Saved speed-up and targeted slow-down experiments to " * OUTPUT_DIR)
end

main()
