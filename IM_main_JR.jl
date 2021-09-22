# Benchmark parameters and functions

using LinearAlgebra, Plots, NLsolve, Interpolations, Optim, Statistics, Distributions

# Parameters
α=2.0
β=1/1.02
γ=0.5
M_min=-2.0
M_max=2.0
z_l=0.9
z_h=1.1
g_l=0.2
g_h=0.25
Pi=[0.25 0.25 0.25 0.25]

# States: ll,lh,hl,hh = 1, 2, 3, 4

# Complete markets solution
path="C:\\Users\\joaor\\Dropbox\\AMSS\\Code\\JR\\"
include(path*"CM_main_JR.jl")

# Incomplete markets functions
include(path*"IM_AMSSFunctions_JR.jl")

# here we construct a 2-d poly approix to the CM V(Z)
#=
A = [-98.15998456008501,-.000001,-.0000001]

function poly_find(A)

    Z_vals = zeros(3)

    Z_vals[1] = Sp_Results[1,1]

        index_0 = findall(x->x==0.0, Sp_Results[:,1])

    Z_vals[2] = Sp_Results[index_0[1],1]

    Z_vals[3] = Sp_Results[end,1]

    V_vals = zeros(3)

    V_vals[1] = Sp_Results[1,11]

    V_vals[2] = Sp_Results[index_0[1],11]

    V_vals[3] = Sp_Results[end,11]

    Err = zeros(3)

    Err = A[1] .+ A[2]*Z_vals + A[3].*Z_vals.^2 - V_vals

    return Err

end


poly_find(A)
results_poly=nlsolve(poly_find,A)
answer_A = results_poly.zero

answer_Amod = answer_A - [.1, 0, .01]

P_fun(x) = answer_Amod[1] + answer_Amod[2]*x + answer_Amod[3]*x^2

Pprim_fun(x) = answer_Amod[2] + 2*answer_Amod[3]*x
=#
#plot(Sp_Results[:,1],[Sp_Results[:,11], P_fun])
#xlabel!("Z")
#ylabel!("Payoff")
#title!("Payoff(Z) values")

#Orig_Z_grid=range(Sp_Results[1,1],Sp_Results[end,1]+0.2,length=length(Sp_Results[:,1])+200)
#V0 = P_fun.(Orig_Z_grid)

Orig_Z_grid=Sp_Results[1:end,1]
V0 = Sp_Results[1:end,11]

interp_V0=LinearInterpolation(Orig_Z_grid,V0;extrapolation_bc=Line())

n_Z=5000
New_Z_grid=range(Orig_Z_grid[1],Orig_Z_grid[end]+0.5,length=n_Z)

V0old=interp_V0(New_Z_grid)


Vf,polf,failedpoints=VfunctionIteration(New_Z_grid)

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
