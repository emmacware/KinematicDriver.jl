## Use Droplets Coagulation for kinematic KiD_driver


function opsplitting_condition(u, t, integrator)
    t % 1.0 == 0.0  # Trigger every 1 second
end

function SD_collision_coalescence_callback(integrator)
    Y = integrator.u
    aux = integrator.p
    aux.microph_variables.SD_Vol .= Y.SD_Vol
    FT = eltype(aux.microph_variables.q_tot)

    coagdata = Droplets.coagulation_run{FT}(aux.droplets_params.Ns)

    (;SD_Vol, SD_Mult) = coalescence_timestep!.(Droplets.Serial(),Droplets.KiD(), 
    aux.microph_variables.SD_Mult,aux.microph_variables.SD_Vol,
    coagdata,aux.droplets_params.settings)
    @. aux.microph_variables.SD_Mult = SD_Mult
    @. aux.microph_variables.SD_Vol = SD_Vol

    @. Y.SD_Vol = aux.microph_variables.SD_Vol
end



function SD_stochasticadv_callback(integrator)
    Y = integrator.u
    aux = integrator.p
    # aux.microph_variables.SD_Vol = Y.SD_Vol
    # dY = integrator.du
    # CO.zero_tendencies!(dY)
    # aux.precip_sinks.SD_Vol .= 0.0
    # aux.scratch.tmp_droplets .= 0.0

    If = CC.Operators.InterpolateC2F()
    # velocity = aux.prescribed_velocity.ρw #/ If(aux.thermo_variables.ρ)
    velocity = aux.thermo_variables.ρ ./ (aux.thermo_variables.ρ.+0.1)
    dz = CC.Fields.Δz_field(aux.thermo_variables.ρ)
    oned_motion_probability = velocity #.* aux.TS.dt ./ dz
    successes = succprob.(Distributions.Binomial.(aux.droplets_params.Ns, oned_motion_probability))

    # map(aux -> move_droplets_to_tmp(aux), aux)

    # transport_temp(aux)
    # should this be upwind?? will this work
    # CC.Operators.RightBiasedCTF(aux.precip_sinks.SD_Vol)
    # CC.Operators.RightBiasedCTF(aux.scratch.tmp_droplets)

    # CC.Operators.RightBiasedFTC(aux.precip_sinks.SD_Vol)
    # CC.Operators.RightBiasedFTC(aux.scratch.tmp_droplets)

    # map(aux -> move_tmp_to_droplets(aux), aux)


    Y.SD_Vol = aux.microph_variables.SD_Vol
end

function move_droplets_to_tmp(aux)

    healthy = findfirst(iszero, aux.microph_variables.SD_Vol) 
    if healthy isa nothing
        healthy = aux.droplets_params.Ns
    end 
    idx = randdperm(healthy)
    aux.microph_variables.SD_Vol[1:healthy] .= aux.microph_variables.SD_Vol[idx]
    aux.microph_variables.SD_Mult[1:healthy] .= aux.microph_variables.SD_Mult[idx]

    transportidx = healthy-aux.scratch.tmp1:healthy
    aux.precip_sinks.SD_Vol[1:aux.scratch.tmp1] .= aux.microph_variables.SD_Vol[transportidx]
    aux.scratch.tmp_droplets[1:aux.scratch.tmp1] .= aux.microph_variables.SD_Mult[transportidx]

    aux.microph_variables.SD_Vol[transportidx] .= 0.0
    aux.microph_variables.SD_Mult[transportidx] .= 0.0

    aux.tmp2 = findfirst(iszero, aux.microph_variables.SD_Vol) 
end

function move_tmp_to_droplets(aux)
    transported = findfirst(iszero, aux.microph_variables.SD_Vol) - 1
    new_idx = aux.scratch.tmp2:aux.scratch.tmp2 + transported
    aux.microph_variables.SD_Vol[new_idx] .= aux.precip_sinks.SD_Vol[1:transported]
    aux.microph_variables.SD_Mult[new_idx] .= aux.scratch.tmp_droplets[1:transported]
end

function transport_temp(aux)
    #shift one cell to the wind direction
    CC.Operators.RightBiasedCTF(aux.microph_variables.SD_Vol)
    CC.Operators.RightBiasedCTF(aux.microph_variables.SD_Mult)
    CC.Operators.RightBiasedFTC(aux.microph_variables.SD_Vol)
    CC.Operators.RightBiasedFTC(aux.microph_variables.SD_Mult)
end