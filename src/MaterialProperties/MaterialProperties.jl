# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

# Maybe not strip the units at some point?
const elementary_charge = Float64(ustrip(u"C", Unitful.q))
const ϵ0 = Float64(ustrip(u"F/m", Unitful.ϵ0))
const kB = Float64(ustrip(u"J/K", Unitful.k))
const me = Float64(ustrip(u"kg", Unitful.me))
const R_gas = Float64(ustrip(u"cal/K/mol", Unitful.R))

const material_properties = Dict{Symbol, NamedTuple}()

abstract type AbstractDriftMaterial end

material_properties[:Vacuum] = (
    E_ionisation = 0.0u"eV",
    f_fano = 0.0,
    ϵ_r = 1.0,
    ρ = 0u"g*cm^-3",
    name = "Vacuum",
    ml = 1.0,
    mt = 1.0
)

abstract type HPGe <: AbstractDriftMaterial end 
Base.Symbol(::Type{HPGe}) = :HPGe
material_properties[:HPGe] = (
    E_ionisation = 2.95u"eV",
    f_fano = 0.129, # https://doi.org/10.1103/PhysRev.163.238
    ϵ_r = 16.0,
    ρ = 5.323u"g*cm^-3",
    name = "High Purity Germanium",
    ml = 1.64,
    mt = 0.0819,
    De = 239u"cm^2/s",
    Dh = 279u"cm^2/s"
)

# These values might just be approximations
abstract type Si <: AbstractDriftMaterial end
Base.Symbol(::Type{Si}) = :Si
material_properties[:Si] = (
    E_ionisation = 3.62u"eV",
    f_fano = 0.11,
    ϵ_r = 11.7,
    ρ = 2.3290u"g*cm^-3",
    name = "Silicon",
    mt = 0.98,
    ml = 0.19,
    De = 39u"cm^2/s",
    Dh = 12u"cm^2/s"
)

abstract type CsPbBr3 <: AbstractDriftMaterial end
Base.Symbol(::Type{CsPbBr3}) = :CsPbBr3
material_properties[:CsPbBr3] = (
    name = "CsPbBr3",
    f_fano = 0.1, # https://doi.org/10.1038/s41566-020-00727-1
    E_ionisation = 6.9u"eV",
    ϵ_r = 16.46,
    ρ = 4.73u"g*cm^-3"
)

material_properties[:Al] = (
    name = "Aluminium",
    ϵ_r = 10.8, # Aluminium Foil
    ρ = 2.6989u"g*cm^-3"
)

material_properties[:Pt] = (
    name = "Platinum",
      ϵ_r = 1,
    ρ = 21.45u"g*cm^-3"
)

material_properties[:LAr] = (
    name = "liquid Argon",
    E_ionisation = 26.0u"eV", # https://doi.org/10.1088/0370-1328/85/6/328
    f_fano = 0.107, # https://doi.org/10.1016/0029-554X(76)90292-5
    ϵ_r = 1.505,
    ρ = 1.396u"g*cm^-3"
)

material_properties[:Cu] = (
    name = "Copper",
    ϵ_r = 20,
    ρ = 8.96u"g*cm^-3"
)

material_properties[:PEN] = (
    name = "Polyethylene Naphthalate",
    ϵ_r = 3.0, # https://topas.com/low-dielectric-constant-plastic-materials-low-permittivity-plastics-topas
    ρ = 1.36u"g*cm^-3"
)

material_properties[:PTFE] = (
    name = "Polytetrafluorethylen",
    ϵ_r = 2.02, # https://topas.com/low-dielectric-constant-plastic-materials-low-permittivity-plastics-topas
    ρ = 2.2u"g*cm^-3"
)

material_properties[:CdZnTe] = (
    name = "Cadmium zinc telluride",
    E_ionisation = 4.64u"eV",
    f_fano = 0.089, # https://doi.org/10.1557/PROC-487-101
    ϵ_r = 10.9,
    ρ = 5.78u"g*cm^-3",
    lifetime_electrons = 3000 * u"ns",             #300 K
    lifetime_holes = 1000 * u"ns",                 #300 K
    De = 25.9 * u"cm^2/s",                         #300 K
    Dh = 1.29 * u"cm^2/s",                         #300 K
    μ_e = 1000 * u"cm^2/(V*s)",
    μ_h = 50 * u"cm^2/(V*s)"
)

material_properties[:CdTe] = (
    name = "Cadmium telluride",
    E_ionisation = 4.42*1e-3 * u"keV",
    f_fano = 0.15,
    ϵ_r = 10.3,
    ρ = 5.85 * u"g*cm^-3",
    lifetime_electrons = 3000.0 * u"ns",
    lifetime_holes = 1300.0 * u"ns",
    De = 27.144593574867194 * u"cm^2/s",
    Dh = 2.326679449274331 * u"cm^2/s",
    μ_e = 1050.0 * u"cm^2/(V*s)",
    μ_h = 90.0 * u"cm^2/(V*s)"
)

#-------H. Allaire--------

material_properties[:Pb] = (
    name = "Lead",
    ϵ_r = 1e6, # high value for lead as for conductor. Usage of lead is shielding.
    ρ = 11.35u"g*cm^-3"
)

# Add new materials above this line
# and just put different spellings into the dict `materials` below

materials = Dict{String, Symbol}(
    "HPGe" => :HPGe,
    "vacuum" => :Vacuum,
    "Vacuum" => :Vacuum,
    "Copper" => :Cu,
    "copper" => :Cu,
    "Al"  => :Al,
    "LAr" => :LAr,
    "CZT" => :CdZnTe,
    "CdTe" => :CdTe,
    "Si" => :Si,
    "Lead" => :Pb
)


for key in keys(material_properties)
    if !haskey(materials, string(key))
        push!(materials, string(key) => key)
    end
end