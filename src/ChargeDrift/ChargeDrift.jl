struct DriftPath{T <: SSDFloat}
    path::Vector{CartesianPoint{T}}
    timestamps::Vector{T}
end

struct EHDriftPath{T <: SSDFloat}
    e_path::Vector{CartesianPoint{T}}
    h_path::Vector{CartesianPoint{T}}
    timestamps_e::Vector{T}
    timestamps_h::Vector{T}
end

function _common_time(dp::EHDriftPath{T})::T where {T <: SSDFloat}
    max(last(dp.timestamps_e), last(dp.timestamps_h))
end
_common_time(dps::Vector{<:EHDriftPath}) =
maximum(_common_time.(dps))

function _common_timestamps(dp::Union{<:EHDriftPath{T}, Vector{<:EHDriftPath{T}}}, Δt) where {T}
    range(zero(Δt), step = Δt, stop = typeof(Δt)(_common_time(dp)) + Δt)
end

@inline function get_velocity_vector(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(velocity_field(pt.x, pt.y, pt.z))
end

@inline function get_velocity_vector(velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, pt::CylindricalPoint{T}) where {T <: SSDFloat}
    return CartesianVector{T}(velocity_field(pt.r, pt.φ, pt.z))
end


function _drift_charges(det::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                        starting_points::VectorOfArrays{CartesianPoint{T}}, energies::VectorOfArrays{T},
                        electric_field::Interpolations.Extrapolation{<:SVector{3}, 3},
                        Δt::RQ; max_nsteps::Int = 2000, diffusion::Bool = false, trapping::Bool = false,
                        σ_init::T = 1.0e-6, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat, RQ <: RealQuantity}

    drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(flatview(starting_points)))
    dt::T = T(to_internal_units(Δt))
    
    drift_path_counter::Int = 0
    
    for (i, start_points) in enumerate(starting_points)
        
        n_hits::Int = length(start_points)
        
        drift_path_e::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
        drift_path_h::Array{CartesianPoint{T}, 2} = Array{CartesianPoint{T}, 2}(undef, n_hits, max_nsteps)
        timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
        timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)
        
        n_e::Int = _drift_charge!( drift_path_e, timestamps_e, det, point_types, start_points, dt, electric_field, Electron; diffusion, trapping, σ_init, verbose,)
        n_h::Int = _drift_charge!( drift_path_h, timestamps_h, det, point_types, start_points, dt, electric_field, Hole; diffusion, trapping, σ_init, verbose,)
    
        for i in eachindex(start_points)
            drift_paths[drift_path_counter + i] = EHDriftPath{T}( drift_path_e[i,1:n_e], drift_path_h[i,1:n_h], timestamps_e[1:n_e], timestamps_h[1:n_h] )
        end
        
        drift_path_counter += n_hits
    end
    
    return drift_paths
end

function modulate_surface_drift(p::CartesianVector{T})::CartesianVector{T} where {T <: SSDFloat}
    return p
end

function modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, vdv::Vector{<:AbstractVirtualVolume{T}})::CartesianVector{T} where {T <: SSDFloat}
    for i in eachindex(vdv)
        if in(pt, vdv[i])
            return modulate_driftvector(sv, pt, vdv[i])
        end
    end
    return sv
end
modulate_driftvector(sv::CartesianVector{T}, pt::CartesianPoint{T}, vdv::Missing) where {T} = sv

@inline function _is_next_point_in_det(pt::AbstractCoordinatePoint{T}, det::SolidStateDetector{T}, point_types::PointTypes{T, 3, S})::Bool where {T <: SSDFloat, S}
    _convert_point(pt, S) in point_types || (pt in det.semiconductor && !(pt in det.contacts))
end

function project_to_plane(v⃗::AbstractArray, n⃗::AbstractArray) #Vector to be projected, #normal vector of plane
    # solve (v⃗+λ*n⃗) ⋅ n⃗ = 0
    # n⃗ = n⃗ ./ norm(n⃗)
    λ = -1 * dot(v⃗, n⃗) / dot(n⃗, n⃗)
    SVector{3,eltype(v⃗)}(v⃗[1] + λ * n⃗[1], v⃗[2] + λ * n⃗[2], v⃗[3] + λ * n⃗[3])
end

function _combine_all_drift!(step_vectors::Vector{CartesianVector{T}}, ∅_vector::CartesianVector{T}, current_pos::Vector{CartesianPoint{T}}, not_done::Vector{Bool},
    drift_path::Array{CartesianPoint{T},2}, electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T}, ::Type{S}, Δσ::T,  
    Δt::T, drift_model::AbstractChargeDriftModel{T}, distribution::Normal{Float64}, istep::Int, point_types::PointTypes{T, 3, S}, τ_p::T, ::Type{CC}, vdv)::Nothing where {T, S, CC <: ChargeCarrier}

    #Loop over indices, with parallel CPU; If particle is still part of the drift
    @inbounds Threads.@threads for n in eachindex(step_vectors)
        
        #_set_to_zero_vector
        step_vectors[n] = ∅_vector

        #_trap_charges!
        if rand(T) > τ_p 
            not_done[n] = false 
        end
        
        #For still moving charges
        if not_done[n] 

            #_add_fieldvector_diffusion!, _add_fieldvector_drift! and _get_driftvectors!
            step_vectors[n] = (getV(SVector{3,T}(step_vectors[n] + get_velocity_vector(electric_field, _convert_point(current_pos[n], S))), drift_model, CC, current_pos[n]) * Δt) + (Δσ * CartesianVector{T}(rand(distribution, 3)))

            #modulate_driftvector
            step_vectors[n] = modulate_driftvector(step_vectors[n], current_pos[n], vdv)

        end

        #_check_and_update_position!
        current_pos[n] += step_vectors[n]

        #If now not normal, crossing a surface
        if (not_done[n]) && !(_is_next_point_in_det(current_pos[n], det, point_types))
            not_done[n] = false
            current_pos[n] = get_crossing_pos(det, drift_path[n,istep-1], current_pos[n])
        end

        #Update drift path
        drift_path[n,istep] = current_pos[n]
    end      
    nothing
end

function _modulate_driftvectors!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, vdv::Vector{V})::Nothing where {T <: SSDFloat, V <: AbstractVirtualVolume{T}}
    for n in eachindex(step_vectors)
        step_vectors[n] = modulate_driftvector(step_vectors[n], current_pos[n], vdv)
    end
    nothing
end
_modulate_driftvectors!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, ::Missing) where {T <: SSDFloat} = nothing

#From Mathieu Benoit, L.A. Hamel, 2009
function _BH_diffusion!(σ_current::T, not_done::Vector{Bool}, α::T, β::T, δ::T)::Tuple{T, T} where {T <: SSDFloat}
    Δσ::T = α * sqrt(2 * (δ + (β * sum(not_done) / σ_current)))
    σ_current = sqrt(σ_current^2 + Δσ^2)
    return Δσ, σ_current 
end

function _drift_charge!(    
                            drift_path::Array{CartesianPoint{T},2},
                            timestamps::Vector{T},
                            det::SolidStateDetector{T},
                            point_types::PointTypes{T, 3, S},
                            startpos::AbstractVector{CartesianPoint{T}},
                            Δt::T,
                            electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3},
                            ::Type{CC};
                            diffusion::Bool = false,
                            trapping::Bool = false,
                            σ_init::T = 1.0e-6,
                            verbose::Bool = true
                        )::Int where {T <: SSDFloat, S, CC <: ChargeCarrier}
        
    #General values
    σ_current::T = σ_init
    n_hits::Int, max_nsteps::Int = size(drift_path)
    ∅_vector::CartesianVector{T} = CartesianVector{T}(0,0,0) 
    distribution::Normal{Float64} = Normal()   

    #Initialize
    drift_path[:,1] = startpos
    timestamps[1] = zero(T)
    last_step::Int = 1

    #Data format preparation
    current_pos::Vector{CartesianPoint{T}} = deepcopy(startpos)
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    not_done::Vector{Bool} = broadcast(pt -> _is_next_point_in_det(pt, det, point_types), startpos)

    #Drift coefficients
    μ_cc::T =
        if CC == Electron
            _parse_value(T, det.semiconductor.material.μ_e, u"cm^2/(V*s)") #Electron mobility   
        else
            _parse_value(T, det.semiconductor.material.μ_h, u"cm^2/(V*s)") #Hole mobility
        end
    drift_model::AbstractChargeDriftModel{T} = det.semiconductor.charge_drift_model    
    verbose && @info "Simulation with a $(CC) mobility coefficient of $(μ_cc)cm^2/(V*s)."

    #Diffusion preparation    
    δ::T = if diffusion #Diffusion coefficient
        if CC == Electron && haskey(det.semiconductor.material, :De)
            _parse_value(T, det.semiconductor.material.De, u"m^2/s")
        elseif CC == Hole && haskey(det.semiconductor.material, :Dh)
            _parse_value(T, det.semiconductor.material.Dh, u"m^2/s")
        else
            zero(T)
        end
    else
        zero(T)
    end

    #Intermediate diffusion coefficients extracted from Mathieu Benoit, L.A. Hamel, 2009
    α::T = if diffusion 
        sqrt(2 * Δt) else zero(T) end  
    β::T = if diffusion 
        (μ_cc * elementary_charge) / (24 * (pi^(3/2)) * ϵ0 * T(det.semiconductor.material.ϵ_r)) else zero(T) end
    verbose && diffusion && @info "Simulation with a $(CC) diffusion coefficient of $(δ)m^2/s."
    
    #Trapping preparation
    τ_p::T = if trapping #Trapping coefficient
        if CC == Electron && haskey(det.semiconductor.material, :lifetime_electrons)
            exp(-Δt/_parse_value(T, det.semiconductor.material.lifetime_electrons, u"s"))
        elseif CC == Hole && haskey(det.semiconductor.material, :lifetime_holes)
            exp(-Δt/_parse_value(T, det.semiconductor.material.lifetime_holes, u"s"))
        else
            one(T)
        end
    else
        one(T)
    end
    verbose && trapping && @info "Simulation with a $(CC) lifetime of $(τ)s."

    #Virtual drift module
    vdv = det.virtual_drift_volumes
    
    #Loop until all charges have drifted or until reaches the maximum number of steps    
    @inbounds for istep in 2:max_nsteps

        #Increment the index number
        last_step += 1

        #Reprocess diffusion coefficient
        if diffusion Δσ, σ_current = _BH_diffusion!(σ_current, not_done, α, β, δ) else Δσ = zero(T) end

        #Process charge drift, including diffusion and trapping
        _combine_all_drift!(step_vectors, ∅_vector, current_pos, not_done, drift_path, electric_field, det, S, Δσ, Δt,
        drift_model, distribution, istep, point_types, τ_p, CC, vdv)

        #Update time
        timestamps[istep] = timestamps[istep-1] + Δt

        #Stop the computation
        if !any(not_done) break end

    end

    return last_step
end

# Point types for charge drift
const CD_ELECTRODE = 0x00
const CD_OUTSIDE = 0x01
const CD_BULK = 0x02
const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

function get_crossing_pos(  det::SolidStateDetector{T}, pt_in::CartesianPoint{T}, pt_out::CartesianPoint{T})::CartesianPoint{T} where {T <: SSDFloat}
    
    # check if the points are already in contacts                    
    if pt_in in det.contacts return pt_in end  
    
    direction::CartesianVector{T} = normalize(pt_out - pt_in)
    crossing_pos::CartesianPoint{T} = pt_out # need undef version for this
    
    # define a Line between pt_in and pt_out
    line::ConstructiveSolidGeometry.Line{T} = ConstructiveSolidGeometry.Line{T}(pt_in, direction)
    tol::T = 5000 * ConstructiveSolidGeometry.csg_default_tol(T)

    # check if the Line intersects with a surface of the Semiconductor
    for surf in ConstructiveSolidGeometry.surfaces(det.semiconductor.geometry)
        for pt in ConstructiveSolidGeometry.intersection(surf, line)
            normal = normalize(ConstructiveSolidGeometry.normal(surf, pt))
            if pt + tol * normal in det.semiconductor &&         # point "before" crossing_pos should be in
                !(pt - tol * normal in det.semiconductor) &&      # point "after" crossing_pos should be out
                0 ≤ (pt - pt_in) ⋅ direction ≤ 1 &&              # pt within pt_in and pt_out
                norm(pt - pt_in) < norm(crossing_pos - pt_in) # pt closer to pt_in that previous crossing_pos
                crossing_pos = pt
            end
        end 
    end

    crossing_pos
end
