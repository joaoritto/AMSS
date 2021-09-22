# To run in a different computer, change the path. Important to use
# double bars in defining the path.
path="C:\\Users\\joaor\\Dropbox\\AMSS\\Code\\JR\\"


include(path*"CM_AMSSFunctions_JR.jl")

using LinearAlgebra, Plots, NLsolve, Interpolations

# Parameters
α=2.0
β=1/1.02
γ=0.5
M_min=-2.0
M_max= 2.0
z_l=0.9
z_h=1.1
g_l=0.2
g_h=0.25
Pi=[0.25 0.25 0.25 0.25]

prelim_sol_Z0,Big_Res = Find_zero()

Guess_sol_Z0 = prelim_sol_Z0[1:5]

ERs_nb(Guess_sol_Z0,0)

function F(x)   # function handle for fsolve
    return ERs_nb(x,0)
end

opt_Z0=nlsolve(F,Guess_sol_Z0)
opt_Lλ_Z0 = opt_Z0.zero

solution_Z0 = zeros(1,11)  # Z=0, lambda, 4 labors, 4 bonds

for count in 1:5
    solution_Z0[1,count+1] = opt_Lλ_Z0[count]
end
for s in 1:4
    solution_Z0[1,s+6] = Bond_LZcm(opt_Lλ_Z0[s+1],0,s)
end
solution_Z0[11] = Pay_Lcm(opt_Lλ_Z0[2:5])


# construct all of the results up until the debt bond binds

ZZ_positive = 0.001:.005:3  # grid counting up

solution_Zpositive = Z_count(ZZ_positive,M_max,opt_Lλ_Z0)

ZZ_negative = -0.001:-.005:-3 # grid counting down

solution_Znegative = Z_count(ZZ_negative,M_min,opt_Lλ_Z0)

Sp_Results = [reverse(solution_Znegative,dims=1); solution_Z0; solution_Zpositive]


p1 = plot(Sp_Results[:,1],Sp_Results[:,3:6])
xlabel!("Z")
ylabel!("Labor")
title!("Labor(Z,s) values")

p2 = plot(Sp_Results[:,1],Sp_Results[:,7:10])
xlabel!("Z")
ylabel!("Bond")
title!("Bond(Z,s) values")

p3 = plot(Sp_Results[:,1],Sp_Results[:,2])
xlabel!("Z")
ylabel!("Lambda")
title!("Lambda(Z) values")

p4 = plot(Sp_Results[:,1],Sp_Results[:,11])
xlabel!("Z")
ylabel!("Payoff")
title!("Payoff(Z) values")

myplot = plot(p1, p2, p3, p4, layout = (2, 2), legend = false)

# Construct guess values and Z at high end

Z_high = Sp_Results[end,1]

Guess_sol_Z_high = zeros(6)
Guess_sol_Z_high[1:5] = Sp_Results[end,2:6]
Guess_sol_Z_high[6] = Z_high

binding_state_high=findmax(Sp_Results[end,7:10])[2]

function H(x)  # function handle for fsolve
    return ERs_bbh(x,Z_high,binding_state_high)
end

Guess_sol_Z_high[4] += - .0015
Guess_sol_Z_high[6] += - .0015

opt_Z_high=nlsolve(H,Guess_sol_Z_high)
opt_Lλ_Z_high = opt_Z_high.zero


# Constructing bottom guess values and Z

Z_low = Sp_Results[1,1]

Guess_sol_Z_low = zeros(6)
Guess_sol_Z_low[1:5] = Sp_Results[1,2:6]
Guess_sol_Z_low[6] = Z_low

binding_state_low=findmin(Sp_Results[1,7:10])[2]

function L(x)  # function handle for fsolve
    return ERs_bbl(x,Z_low,binding_state_low)
end

Guess_sol_Z_low[5] += - .001
Guess_sol_Z_low[6] += + .001

opt_Z_low=nlsolve(L,Guess_sol_Z_low)
opt_Lλ_Z_low = opt_Z_low.zero

# Add to table Sp_Results, of Z in which gov debt binds
Sp_Results[1,3:6]=opt_Lλ_Z_low[2:5]
Sp_Results[1,11]=Pay_Lcm(opt_Lλ_Z_low[2:5])
Sp_Results[end,3:6]=opt_Lλ_Z_high[2:5]
Sp_Results[end,11]=Pay_Lcm(opt_Lλ_Z_high[2:5])
