## 2020-12-11
## Jani V Anttila

using JSON;

## Flux struct - used for bookkeeping for output purposes

mutable struct Flux
    fluxin::Array{Float64,1}
    fluxside::Array{Float64,1}
    fluxtop_d::Array{Float64,1}
    fluxtop_a::Array{Float64,1}
end

function make_flux(in_ny::Integer,in_nr::Integer)::Flux
    fin = zeros(in_nr)
    fs = zeros(in_ny)
    ftd = zeros(in_nr)
    fta = zeros(in_nr)
    f = Flux(fin,fs,ftd,fta)
    return f
end

function flux_totals(f::Flux)::Array{Float64,1}
    res = zeros(4)
    res[1] = sum(f.fluxin)
    res[2] = sum(f.fluxside)
    res[3] = sum(f.fluxtop_d)
    res[4] = sum(f.fluxtop_a)
    return res
end

Base.show(io::IO, f::Flux) = print(io, "in    ",(flux_totals(f))[1],
                                   " \nside  ",(flux_totals(f))[2],
                                   " \ntop_d ",(flux_totals(f))[3],
                                   " \ntop_a ",(flux_totals(f))[4])

## state struct - holds the state variables of the tree

mutable struct State
    cair::Array{Float64,2} ## conc. in gaseous phase
    cwat::Array{Float64,2} ## conc. in aqueous phase
    nair::Array{Float64,2} ## n [mol] in gaseous phase
    nwat::Array{Float64,2} ## n [mol] in aqueous phase
    time::Float64
    ny::Integer
    nr::Integer
end

function make_state(in_ny::Integer,in_nr::Integer)::State
    ca = zeros(in_ny,in_nr)
    cw = zeros(in_ny,in_nr)
    na = zeros(in_ny,in_nr)
    nw = zeros(in_ny,in_nr)
    s = State(ca,cw,na,nw,0.0,in_ny,in_nr)
end

function total_storage(s::State)
    nsum = 0.0
    for i = 1:s.ny
        for j = 1:s.nr
            nsum += s.nair[i,j]
            nsum += s.nwat[i,j]
        end
    end
    return nsum
end

Base.show(io::IO, s::State) = print(io, "time ",s.time,
                                    " ny ", s.ny, " nr ", s.nr,
                                    " total ", total_storage(s))

## derivative state - holds the state change

mutable struct Derivative
    qair::Array{Float64,2}    ## radial diffusion
    qadv::Array{Float64,2}    ## advection
    axdiff::Array{Float64,2}  ## axial diffusion
    towater::Array{Float64,2} ## phase equilibration
end

function make_derivative(in_ny::Integer,in_nr::Integer)::Derivative
    qa = zeros(in_ny,in_nr)
    qd = zeros(in_ny,in_nr)
    axd = zeros(in_ny,in_nr)
    tow = zeros(in_ny,in_nr)
    d = Derivative(qa,qd,axd,tow)
end

mutable struct Cumulant
    fluxin::Float64
    fluxside::Float64
    fluxtop_d::Float64
    fluxtop_a::Float64
end

Base.show(io::IO, c::Cumulant) = print(io, "in    ",c.fluxin,
                                       " \nside  ",c.fluxside,
                                       " \ntop_d ",c.fluxtop_d,
                                       " \ntop_a ",c.fluxtop_a)

mutable struct Tree
    volume::Array{Float64,2}
    volumeair::Array{Float64,2}
    volumewat::Array{Float64,2}
    crossect::Array{Float64,2}
    ker1::Array{Float64,2}
    ker2::Array{Float64,2}
    ny::Integer
    nr::Integer
    initialised::Bool
end

function make_tree(in_ny::Integer,in_nr::Integer)::Tree
    vol = zeros(in_ny,in_nr)
    vair = zeros(in_ny,in_nr)
    vwat = zeros(in_ny,in_nr)
    csect = zeros(in_ny,in_nr)
    k1 = zeros(in_ny,in_nr)
    k2 = zeros(in_ny,in_nr)
    t = Tree(vol,vair,vwat,csect,k1,k2,in_ny,in_nr,false)
end

function crossect_area(idx::Integer,dr::Float64)::Float64
    if idx == 1
        return pi*dr*dr
    else
        return pi*(idx * dr)^2 - pi*((idx - 1) * dr)^2
    end
end

function initialise_tree_fix_dr(t::Tree,p::Dict{String,Any})
    for i = 1:t.ny
        for j = 1:t.nr
            t.crossect[i,j] = crossect_area(j,p["dr"])
        end
    end
    for i = 1:t.ny
        t.volume[i,1] = p["dy"] * t.crossect[i,1]
        t.volumeair[i,1] = p["fractair"] * t.volume[i,1]
        t.volumewat[i,1] = p["fractwat"] * t.volume[i,1]
        t.ker1[i,1] = Inf
        t.ker2[i,1] = log(2.0 * p["dr"] / p["dr"])
        for j = 2:t.nr
            t.volume[i,j] = p["dy"] * t.crossect[i,j]
            t.volumeair[i,j] = p["fractair"] * t.volume[i,j]
            t.volumewat[i,j] = p["fractwat"] * t.volume[i,j]
            t.ker1[i,j] = log((p["dr"] * j) / (p["dr"] * (j-1)))
            t.ker2[i,j] = log((p["dr"] * (j+1)) / (p["dr"] * j))
        end
    end
    t.initialised = true;
    return nothing;
end

function radial_diffusion(t::Tree,s::State,d::Derivative,f::Flux,
                          p::Dict{String,Any})
    ## set derivative elements
    for i = 1:t.ny
        for j = 1:t.nr
	    if j == 1
	        ## innermost radial element
	        d.qair[i,j] = (-2.0 * pi * p["diff_r"] * (s.cair[i,j] - s.cair[i,j+1])
                               / t.ker2[i,j])
	    elseif j == t.nr
	        ## outermost radial element
	        d.qair[i,j] = ((2.0 * pi * p["diff_r"] * (s.cair[i,j-1] - s.cair[i,j])
                                / t.ker1[i,j]) -
			       (2.0 * pi * p["diff_b"] * (s.cair[i,j] - p["camb"])
                                / t.ker2[i,j]))
	    else
	        d.qair[i,j] = ((2.0 * pi * p["diff_r"] * (s.cair[i,j-1] - s.cair[i,j])
                                / t.ker1[i,j]) -
			       (2.0 * pi * p["diff_r"] * (s.cair[i,j] - s.cair[i,j+1])
                                / t.ker2[i,j]))
	    end
        end
    end
    ## set flux elements
    for i = 1:t.ny
	f.fluxside[i] = (2.0 * pi * p["diff_b"] * (s.cair[i,t.nr] - p["camb"])
                         / t.ker2[i,t.nr]);
    end
    return nothing;
end

function axial_advection(t::Tree,s::State,d::Derivative,f::Flux,
                         p::Dict{String,Any})
    ## set derivative elements
    for i = 1:t.ny
        for j = 1:t.nr
	    if j <= p["nrs"]
	        ## zero inside the duramen
	        d.qadv[i,j] = 0.0;
	    else
	        if i == 1
	            ## bottom element layer receives CH4 in water
	            d.qadv[i,j] = ((p["vel"] * p["fractwat"] * p["csoil"] * t.crossect[i,j]) -
			           (p["vel"] * s.cwat[i,j] * t.crossect[i,j]));
	            f.fluxin[j] = (p["vel"] * p["fractwat"] * p["csoil"] * t.crossect[i,j]);
	        elseif i == t.ny
	            ## the topmost elements simply lose CH4 in water
	            d.qadv[i,j] = ((p["vel"] * s.cwat[i-1,j] * t.crossect[i-1,j]) -
			           (p["vel"] * s.cwat[i,j] * t.crossect[i,j]));
	        else
	            d.qadv[i,j] = ((p["vel"] * s.cwat[i-1,j] * t.crossect[i-1,j]) -
			           (p["vel"] * s.cwat[i,j] * t.crossect[i,j]));
	        end
	    end
        end
    end
    ## set flux elements
    for j = 1:t.nr
        if j <= p["nrs"]
	    f.fluxtop_a[j] = 0.0;
        else
	    f.fluxtop_a[j] = p["vel"] * s.cwat[t.ny,j] * t.crossect[t.ny,j];
        end
    end
    return nothing;
end

function axial_diffusion(t::Tree,s::State,d::Derivative,f::Flux,
                         p::Dict{String,Any})
    ## set derivative elements
    for i = 1:t.ny
        for j = 1:t.nr
	    if i == 1
	        ## without axial diffusion in
	        d.axdiff[i,j] = -(p["diff_a"] * t.crossect[i,j] *
			          (s.cair[i,j] - s.cair[i+1,j]) / p["dy"]);
	        ## with axial diffusion in
	        ## axdiff[i,j]=(p["diff_a"]*crossect[j]*(p["csoil"]-s.cair[i,j]) / p["dy"]) -
	        ##         (p["diff_a"]*crossect[j]*(s.cair[i,j] - s.cair[i+1,j]) / p["dy"]);
	    elseif i == t.ny
	        d.axdiff[i,j] = ((p["diff_a"] * t.crossect[i,j] *
			          (s.cair[i-1,j] - s.cair[i,j]) / p["dy"]) -
			         (p["diff_a"] * t.crossect[i,j] *
			          (s.cair[i,j] - p["camb"]) / p["dy"]));
	    else
	        d.axdiff[i,j] = ((p["diff_a"] * t.crossect[i,j] *
			          (s.cair[i-1,j] - s.cair[i,j]) / p["dy"]) -
			         (p["diff_a"] * t.crossect[i,j] *
			          (s.cair[i,j] - s.cair[i+1,j]) / p["dy"]));
	    end
        end
    end
    ## set flux elements
    for j = 1:t.nr
    	f.fluxtop_d[j] = (p["diff_a"] * t.crossect[t.ny,j] *
                          (s.cair[t.ny,j] - p["camb"]) / p["dy"]);
    end
    return nothing;
end

function phase_equilibration(t::Tree,s::State,d::Derivative,p::Dict{String,Any})
    for i = 1:t.ny
        for j = 1:t.nr
	    d.towater[i,j] = ((p["henry"] * s.cair[i,j] - s.cwat[i,j]) *
                              t.volumewat[i,j] * p["gamma"]);
        end
    end
    return nothing;
end

function update_cumulants(c::Cumulant,a::Array{Float64,1},dt::Float64)
    c.fluxin    += dt * a[1];
    c.fluxside  += dt * a[2];
    c.fluxtop_d += dt * a[3];
    c.fluxtop_a += dt * a[4];
    return nothing;
end

function set_c_from_n(t::Tree,s::State)
    for i = 1:t.ny
        for j = 1:t.nr
            s.cair[i,j] = s.nair[i,j]/t.volumeair[i,j];
            s.cwat[i,j] = s.nwat[i,j]/t.volumewat[i,j];
        end
    end
    return nothing;
end

## todo : check that tree is initialised
function euler_step(t::Tree,s::State,d::Derivative,f::Flux,
                    c::Cumulant,p::Dict{String,Any})
    for step = 1:p["res"]
	radial_diffusion(t,s,d,f,p);
	axial_advection(t,s,d,f,p);
	axial_diffusion(t,s,d,f,p);
	phase_equilibration(t,s,d,p);
	##update_cumulants(p);
	for i = 1:t.ny
	    for j = 1:t.nr
	        s.nair[i,j] += p["dt"]*(d.qair[i,j] +
					d.axdiff[i,j] -
					d.towater[i,j]);
	        s.nwat[i,j] += p["dt"]*(d.qadv[i,j] +
					d.towater[i,j]);
	        s.cair[i,j] = s.nair[i,j]/t.volumeair[i,j];
	        s.cwat[i,j] = s.nwat[i,j]/t.volumewat[i,j];
	    end
	end
        update_cumulants(c,flux_totals(f),p["dt"]);
        s.time += p["dt"]
    end
    return nothing;
end

