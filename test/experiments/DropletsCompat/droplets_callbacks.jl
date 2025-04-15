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

    # Droplets.static_droplets_attributes{FT,Nsd}(aux.microph_variables.SD_Mult,aux.microph_variables.SD_Vol)

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


    down_to_face = CC.Operators.RightBiasedC2F()
    up_to_face = CC.Operators.LeftBiasedC2F()
    down_to_center = CC.Operators.RightBiasedF2C()
    up_to_center = CC.Operators.LeftBiasedF2C()

    If = CC.Operators.InterpolateC2F()

    dz = CC.Fields.Δz_field(aux.thermo_variables.ρ)
    dt = aux.TS.dt

    uppervel = similar(aux.thermo_variables.ρ)
    lowervel = similar(aux.thermo_variables.ρ)

    @. uppervel = down_to_center.(aux.prescribed_velocity.ρw.components.data.:1 / If.(aux.thermo_variables.ρ))
    @. lowervel = up_to_center.(aux.prescribed_velocity.ρw.components.data.:1 / If.(aux.thermo_variables.ρ))

    @. uppervel *= dt / dz
    @. lowervel *= dt / dz

    (;moveup,movedown) = advection_trials.(aux.droplets_params.Ns, uppervel,lowervel)

    (;SD_Vol,SD_Mult,transport_vol_up,transport_vol_down,transport_mult_up,transport_mult_down
        ) = map(movedropletstotmp, aux.microph_variables.SD_Vol, 
            aux.microph_variables.SD_Mult,moveup, movedown)

    @. transport_vol_up = up_to_center.(up_to_face.(transport_vol_up))
    @. transport_vol_down = down_to_center.(down_to_face.(transport_vol_down))
    @. transport_mult_up = up_to_center.(up_to_face.(transport_mult_up))
    @. transport_mult_down = down_to_center.(down_to_face.(transport_mult_down))


    (;SD_Vol,SD_Mult) = map(finishmovedroplets, SD_Vol,SD_Mult, transport_vol_up,
        transport_vol_down,transport_mult_up,
        transport_mult_down)

    @. Y.SD_Vol = SD_Vol
    @. aux.microph_variables.SD_Vol = SD_Vol
    @. aux.microph_variables.SD_Mult = SD_Mult

end

function advection_trials(Nsd::Int, movedown::FT,moveup::FT) where {FT}
    movedownprob = movedown < 0.0 ? movedown : FT(0.0)
    moveupprob = moveup > 0.0 ? moveup : FT(0.0)
    stayprob = FT(1.0) - movedownprob - moveupprob
    probabilities = SVector{3,FT}(
        movedownprob,
        moveupprob,
        stayprob,
    )
    moveup,movedown,stay = FT.(rand(Distributions.Multinomial(Nsd, probabilities)))
    return (;moveup,movedown)
end

function movedropletstotmp(SD_Vol::SVector{Nsd,FT}, SD_Mult::SVector{Nsd,FT}, transport_up::FT,transport_down) where {FT, Nsd}
    Ns::Int = length(SD_Vol)
    unhealthy::Int = findfirst(iszero, SD_Vol) === nothing ? Ns : findfirst(iszero, SD_Vol)

    # Need to be shuffled but that shouldn't happen here

    transportupidx::UnitRange{Int} = (unhealthy - Int(transport_up)):(unhealthy-1)
    transportdownidx::UnitRange{Int} = (unhealthy -Int(transport_up)-Int(transport_down)):(unhealthy-Int(transport_up)-1)

    transport_vol_up = SVector{Nsd,FT}(i in transportupidx ? SD_Vol[i] : FT(0.0) for i in 1:Ns)
    transport_mult_up = SVector{Nsd,FT}(i in transportupidx ? SD_Mult[i] : FT(0.0) for i in 1:Ns)
    transport_vol_down = SVector{Nsd,FT}(i in transportdownidx ? SD_Vol[i] : FT(0.0) for i in 1:Ns)
    transport_mult_down = SVector{Nsd,FT}(i in transportdownidx ? SD_Mult[i] : FT(0.0) for i in 1:Ns)

    SD_Vol::SVector{Nsd,FT} = SVector{Nsd,FT}(i < (unhealthy-Int(transport_down-transport_up)) ? FT(SD_Vol[i]) : FT(0.0) for i in 1:Ns)
    SD_Mult::SVector{Nsd,FT} = SVector{Nsd,FT}(i < (unhealthy-Int(transport_down-transport_up)) ? FT(SD_Mult[i]) : FT(0.0) for i in 1:Ns)

    return (;SD_Vol=SD_Vol, SD_Mult=SD_Mult, transport_vol_up=transport_vol_up,
        transport_vol_down=transport_vol_down, transport_mult_up=transport_mult_up,
        transport_mult_down=transport_mult_down)
end

function finishmovedroplets(SD_Vol::SVector{Nsd,FT}, SD_Mult::SVector{Nsd,FT}, transport_vol_up::SVector{Nsd,FT},
    transport_vol_down::SVector{Nsd,FT},transport_mult_up::SVector{Nsd,FT},
    transport_mult_down::SVector{Nsd,FT}) where {FT,Nsd}
    Ns::Int = length(SD_Vol)
    first_zero_cell::Int = findfirst(iszero, SD_Vol) === nothing ? Ns : findfirst(iszero, SD_Vol)
    first_zero_tranport_up::Int = findfirst(iszero, transport_vol_up) === nothing ? Ns : findfirst(iszero, transport_vol_up)
    first_zero_tranport_down::Int = findfirst(iszero, transport_vol_down) === nothing ? Ns : findfirst(iszero, transport_vol_down)

    SD_Vol::SVector{Nsd, FT} = SVector{Nsd, FT}([
        i < first_zero_cell ? SD_Vol[i] :
        first_zero_cell < i <= first_zero_cell + first_zero_tranport_up ? transport_vol_up[i - first_zero_cell] :
        first_zero_cell + first_zero_tranport_up < i <= first_zero_cell + first_zero_tranport_up + first_zero_tranport_down ? transport_vol_down[i - first_zero_cell - first_zero_tranport_up] :
        FT(0.0) for i in 1:Nsd
    ])
    SD_Mult::SVector{Nsd, FT} = SVector{Nsd, FT}([
        i < first_zero_cell ? SD_Mult[i] :
        first_zero_cell < i <= first_zero_cell + first_zero_tranport_up ? transport_mult_up[i - first_zero_cell] :
        first_zero_cell + first_zero_tranport_up < i <= first_zero_cell + first_zero_tranport_up + first_zero_tranport_down ? transport_mult_down[i - first_zero_cell - first_zero_tranport_up] :
        FT(0.0) for i in 1:Nsd
    ])

    return (;SD_Vol=SD_Vol,SD_Mult=SD_Mult)
end