include("test.jl")
using JLD, DataFrames
fires = load("/Volumes/Data/Data/STPP/Forest Fire/Data/bridger_teton.jld")["fires"];

include("thomas_pp.jl")
include("utility_functions.jl")

x = fires[:longitude];
x -= minimum(x);
x /= maximum(x);
y = fires[:latitude];
y -= minimum(y);
y /= maximum(y);
x = Vector(x);
y = Vector(y);
id = find(fires[:fire_year].==2004);
out_thomas = thomas_mple(x[id], y[id], 15., 8., .01, process=3, showplot=false);
out_inv = inverse_power_mple(x[id], y[id], 30., 8., 1.5, .005, .5, skipi=100, process=3, showplot=false);
out_genA = thomas_genA_mple(x[id], y[id], 4., 10., .3, .01,.02, .5, skipi=100, process=3, showplot=false);
out_genB = thomas_genB_mple(x[id], y[id], 3., 10., 30., .02, .03, process=3, showplot=false);
out_genC = thomas_genC_mple(x[id], y[id], 40., 10., 15., 3., .02, .03, process=3, showplot=false);

AIC_vals = zeros(5);
BIC_vals = zeros(5);
AIC_vals[1] = out_thomas.AIC;
AIC_vals[2] = out_inv.AIC;
AIC_vals[3] = out_genA.AIC;
AIC_vals[4] = out_genB.AIC;
AIC_vals[5] = out_genC.AIC;

BIC_vals[1] = out_thomas.BIC;
BIC_vals[2] = out_inv.BIC;
BIC_vals[3] = out_genA.BIC;
BIC_vals[4] = out_genB.BIC;
BIC_vals[5] = out_genC.BIC;

x1,y1,xp,yp = thomasPP(out_thomas.pa_sim[end,1],out_thomas.pa_sim[end,2],out_thomas.pa_sim[end,3]);
