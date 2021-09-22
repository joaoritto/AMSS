# Benchmark parameters and functions

addprocs(2-length(procs()))

using Distributed,SharedArrays, Plots, Distributions
@everywhere using LinearAlgebra, NLsolve, Interpolations, Statistics



# Parameters
@everywhere α=2.0
@everywhere β=1/1.02
@everywhere γ=0.5
@everywhere M_min=-2.0
@everywhere M_max=2.0
@everywhere z_l=0.9
@everywhere z_h=1.1
@everywhere g_l=0.2
@everywhere g_h=0.25
@everywhere Pi=[0.25 0.25 0.25 0.25]

# States: ll,lh,hl,hh = 1, 2, 3, 4

# Complete markets solution
path="C:\\Users\\joaor\\Dropbox\\AMSS\\Code\\JR\\"
include(path*"CM_main_JR_parallel.jl")

# Incomplete markets functions
include(path*"IM_AMSSFunctions_JR_parallel.jl")


Orig_Z_grid=Sp_Results[1:end,1]
V0 = Sp_Results[1:end,11]

interp_V0=LinearInterpolation(Orig_Z_grid,V0;extrapolation_bc=Line())

n_Z=5000
New_Z_grid=range(Orig_Z_grid[1],Orig_Z_grid[end]+0.5,length=n_Z)

V0old=interp_V0(New_Z_grid)



Vf,polf,failedpoints=VfunctionIteration(New_Z_grid,Sp_Results,interp_V0)

(λZ,pol_L,pol_B,pol_Z)=polf

p1=plot(New_Z_grid,Vf)
plot!(Sp_Results[:,1],Sp_Results[:,11])
plot!(Orig_Z_grid,V0)
display(p1)



# Simulation

T=1000
debt=zeros(T)
state=zeros(Int64,T)
labor=zeros(T)
consumption=zeros(T)
taxrate=zeros(T)
govexpenditure=zeros(T)
Zi=zeros(Int64,T)
Zi[1]=findmin(New_Z_grid.^2)[2]
for t in 1:T
    s=ceil(rand(Uniform(0,4)))
    zs,gs=ExogStates(s)
    state[t]=s
    debt[t]=pol_B[Zi[t]]
    labor[t]=pol_L[Zi[t],state[t]]
    consumption[t]=zs*labor[t]-gs
    govexpenditure[t]=gs
    taxrate[t]=1-((labor[t]^γ)/(zs*consumption[t]^(-α)))
    if t<T
        Zi[t+1]=findmin((ones(n_Z)*pol_Z[Zi[t],state[t]]-New_Z_grid).^2)[2]
    end
end



T=1000
debtCM=zeros(T)
laborCM=zeros(T)
consumptionCM=zeros(T)
taxrateCM=zeros(T)
govexpenditureCM=zeros(T)
Zi0=findmin(Orig_Z_grid.^2)[2]
for t in 1:T
    zs,gs=ExogStates(state[t])
    debtCM[t]=Sp_Results[Zi0,6+state[t]]
    laborCM[t]=Sp_Results[Zi0,2+state[t]]
    consumptionCM[t]=zs*laborCM[t]-gs
    govexpenditureCM[t]=gs
    taxrateCM[t]=1-((laborCM[t]^γ)/(zs*consumptionCM[t]^(-α)))
end

p2=plot(consumption)
plot!(consumptionCM)
p3=plot(debt)
plot!(debtCM)
p4=plot(taxrate)
plot!(taxrateCM)
p5=plot(labor)
plot!(laborCM)
#p6=plot(Zstate)

finalplot=plot(p2,p3,p4,p5,layout=(2,2))
