## 2021-11-01
## Jani V Anttila

include("stem-gas-transport.jl")

read_par = read("example.json", String)
par = JSON.parse(read_par)
par["dy"] = par["height"]/par["ny"]

t0 = make_tree(par["ny"],par["nr"])
s0 = make_state(par["ny"],par["nr"])
d0 = make_derivative(par["ny"],par["nr"])
f0 = make_flux(par["ny"],par["nr"])
c0 = Cumulant(0.0,0.0,0.0,0.0)

initialise_tree_fix_dr(t0,par)

euler_step(t0,s0,d0,f0,c0,par)
