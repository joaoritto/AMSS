

@everywhere function IM_FOC1(Labor,Z,λold,Z_grid)

    interp_derV=LinearInterpolation(Z_grid,-λold;extrapolation_bc=Line())

    cum1 = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        if cons <= 0

            println("someone did something very bad here")
            error=1000*ones(4,1)
            return error
        end

        cum1 += cons^(-α)*Pi[s]

    end

    B = Z / cum1

    if B<M_max && B>M_min

        nZ = zeros(4)
        μ = zeros(4)
        num_λ=0
        den_λ=0

        for s = 1:4

            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

            μ[s]=interp_derV(nZ[s])

            num_λ+=μ[s]*cons^(-α)*Pi[s]
            den_λ+=cons^(-α)*Pi[s]
        end

        λ=-num_λ/den_λ

        error=zeros(4)
        for s in 1:4
            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            error[s]=cons^(-α)*zs-labor^γ-(μ[s]+λ)*α*cons^(-α-1)*zs*B-μ[s]*
            ((1-α)*cons^(-α)*zs-(1+γ)*labor^γ)

        end


    else
        if B<M_min
            B=M_min
        elseif B>M_max
            B=M_max
        end

        nZ = zeros(4)
        μ = zeros(4)
        error=zeros(4)

        for s = 1:4
            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

            μ[s]=interp_derV(nZ[s])

            if s==1
                num_λ=cons^(-α)*zs-labor^γ-μ[s]*α*cons^(-α-1)*zs*B-μ[s]*((1-α)*cons^(-α)*zs-(1+γ)*labor^γ)
                den_λ=B*α*cons^(-α-1)*zs
                λ=num_λ/den_λ
            else
                error[s]=cons^(-α)*zs-labor^γ-μ[s]*α*cons^(-α-1)*zs*B-μ[s]*((1-α)*cons^(-α)*zs-(1+γ)*labor^γ)-
                λ*B*α*cons^(-α-1)*zs
            end
        end

        for s = 1:4

            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            error[1]+=B*cons^(-α)*Pi[s]
        end
        error[1]=error[1]-Z
    end

    totalerror=sum(abs.(error))

    return error
end

@everywhere function IM_λ(Labor,Z,λold,Z_grid)

    interp_derV=LinearInterpolation(Z_grid,-λold;extrapolation_bc=Line())

    cum1 = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        cum1 += cons^(-α)*Pi[s]

    end

    B = Z / cum1

    if B<M_max && B>M_min

        nZ = zeros(4)
        μ = zeros(4)
        num_λ=0
        den_λ=0

        for s = 1:4

            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

            μ[s]=interp_derV(nZ[s])

            num_λ+=μ[s]*cons^(-α)*Pi[s]
            den_λ+=cons^(-α)*Pi[s]
        end

        λ=-num_λ/den_λ

    else
        if B<M_min
            B=M_min
        elseif B>M_max
            B=M_max
        end

        nZ = zeros(4)
        μ = zeros(4)
        error=zeros(4)

        s=1
        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

        μ[s]=interp_derV(nZ[s])

        num_λ=cons^(-α)*zs-labor^γ-μ[s]*α*cons^(-α-1)*zs*B-μ[s]*((1-α)*cons^(-α)*zs-(1+γ)*labor^γ)
        den_λ=B*α*cons^(-α-1)*zs
        λ=num_λ/den_λ
    end
    return λ
end

@everywhere function IM_B(Labor,Z,λold,Z_grid)

    interp_derV=LinearInterpolation(Z_grid,-λold;extrapolation_bc=Line())

    cum1 = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        cum1 += cons^(-α)*Pi[s]

    end

    B = min(max(Z / cum1,M_min),M_max)
    return B
end

@everywhere function IM_nZ(Labor,Z,λold,Z_grid)

    interp_derV=LinearInterpolation(Z_grid,-λold;extrapolation_bc=Line())

    cum1 = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        cum1 += cons^(-α)*Pi[s]

    end

    B = Z / cum1

    if B<M_max && B>M_min

        nZ = zeros(4)
        μ = zeros(4)
        num_λ=0
        den_λ=0

        for s = 1:4

            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

            μ[s]=interp_derV(nZ[s])

            num_λ+=μ[s]*cons^(-α)*Pi[s]
            den_λ+=cons^(-α)*Pi[s]
        end

        λ=-num_λ/den_λ

    else
        if B<M_min
            B=M_min
        elseif B>M_max
            B=M_max
        end

        nZ = zeros(4)
        μ = zeros(4)
        error=zeros(4)

        for s = 1:4

            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

        end
    end
    return nZ
end

@everywhere function IM_Payoff(Labor,Z,Vold,Z_grid)

    interp_V=LinearInterpolation(Z_grid,Vold;extrapolation_bc=Line())

    cum1 = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        cum1 += cons^(-α)*Pi[s]

    end

    B = Z / cum1
    Payoff=0

    if B<M_max && B>M_min

        nZ = zeros(4)


        for s = 1:4

            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

            VnZ=interp_V(nZ[s])

            Payoff+=(cons^(1-α)/(1-α)-labor^(1+γ)/(1+γ)+β*VnZ)*Pi[s]

        end

    else
        if B<M_min
            B=M_min
        elseif B>M_max
            B=M_max
        end

        nZ = zeros(4)
        μ = zeros(4)
        error=zeros(4)

        for s = 1:4
            zs, gs = ExogStates(s)

            labor = Labor[s]

            cons = zs*labor - gs

            nZ[s] = ( cons^(-α)*B - (cons^(1-α) - labor^(1+γ)) )/β

            VnZ=interp_V(nZ[s])

            Payoff+=(cons^(1-α)/(1-α)-labor^(1+γ)/(1+γ)+β*VnZ)*Pi[s]
        end
    end

    return Payoff
end

function VfunctionIteration(Z_grid,CM_Results,interp_V0)

    @eval @everywhere Z_grid=$Z_grid
    @eval @everywhere CM_Results=$CM_Results
    @eval @everywhere interp_V0=$interp_V0
    @everywhere ssizeZ=size(Z_grid)[1]


    Vold=interp_V0(Z_grid)
    λold=zeros(ssizeZ)

    for j in 2:ssizeZ-1
        if j<ssizeZ
            λold[j]=-(Vold[j+1]-Vold[j-1])/(Z_grid[j+1]-Z_grid[j-1])
        else
            if j==1
                λold[j]=-(Vold[j+1]-Vold[j])/(Z_grid[j+1]-Z_grid[j])
            elseif j==ssizeZ
                λold[j]=-(Vold[j]-Vold[j-1])/(Z_grid[j]-Z_grid[j-1])
            end
        end
    end

    Vnew=SharedArray{Float64}(ssizeZ)
    λnew=SharedArray{Float64}(ssizeZ)
    Labor_sol=SharedArray{Float64}(ssizeZ,4)
    B_sol=SharedArray{Float64}(ssizeZ)
    nZ_sol=SharedArray{Float64}(ssizeZ,4)

    ϵ=1e-10
    error=1000
    iter=0

    while error>ϵ && iter<2000
        iter=iter+1

        if iter==1
            midpoint=ceil(Int64,ssizeZ/2)
            Zi_points=[midpoint:1:ssizeZ;midpoint-1:-1:1]
            failedpoints=[]
            for j in 1:ssizeZ
                Zi=Zi_points[j]
                Z=Z_grid[Zi]
                IM_FOC(Labor)=IM_FOC1(Labor,Z,λold,Z_grid)
                Lab_guess=zeros(4)
                for s in 1:4
                    if j==1
                        lab_interp=LinearInterpolation(Sp_Results[:,1],Sp_Results[:,2+s],extrapolation_bc=Line())
                        Lab_guess[s] = lab_interp(Z)
                    else
                        Lab_guess[s]=Labor_sol[Zi_points[j-1],s]
                    end
                end
                opt = nlsolve(IM_FOC,Lab_guess,iterations=10000)
                if maximum(IM_FOC(opt.zero))>1e-6
                    println("bad zero, Zi=",Zi)
                    push!(failedpoints,Zi)
                #else
                #    println("good zero, Zi=",Zi)
                end

                Labor_sol[Zi,:]=opt.zero
                Vnew[Zi]=IM_Payoff(opt.zero,Z,Vold,Z_grid)
                λnew[Zi]=IM_λ(opt.zero,Z,λold,Z_grid)
                B_sol[Zi]=IM_B(opt.zero,Z,λold,Z_grid)
                nZ_sol[Zi,:]=IM_nZ(opt.zero,Z,λold,Z_grid)
            end
            if length(failedpoints)>1000 && failedpoints[end]>1
                for j in length(failedpoints):-1:1
                    Zi=failedpoints[j]
                    Z=Z_grid[Zi]
                    Lab_guess=Labor_sol[failedpoints[j]-1,:]


                    IM_FOC(Labor)=IM_FOC1(Labor,Z,λold,Z_grid)
                    opt = nlsolve(IM_FOC,Lab_guess,iterations=10000)
                    if maximum(abs.(IM_FOC(opt.zero)))>1e-6
                    else
                        setdiff!(failedpoints,Zi)
                    end
                    Labor_sol[Zi,:]=opt.zero
                    Vnew[Zi]=IM_Payoff(opt.zero,Z,Vold,Z_grid)
                    λnew[Zi]=IM_λ(opt.zero,Z,λold,Z_grid)
                    B_sol[Zi]=IM_B(opt.zero,Z,λold,Z_grid)
                    nZ_sol[Zi,:]=IM_nZ(opt.zero,Z,λold,Z_grid)
                end
            end
        else
            @eval @everywhere λold=$λold
            @eval @everywhere Vold=$Vold

            @eval @everywhere Labor_solold=$Labor_sol

            Vnew=SharedArray{Float64}(ssizeZ)
            λnew=SharedArray{Float64}(ssizeZ)
            Labor_sol=SharedArray{Float64}(ssizeZ,4)
            B_sol=SharedArray{Float64}(ssizeZ)
            nZ_sol=SharedArray{Float64}(ssizeZ,4)

            @sync @distributed for Zi in 1:ssizeZ
                Z=Z_grid[Zi]
                IM_FOC(Labor)=IM_FOC1(Labor,Z,λold,Z_grid)
                Lab_guess=zeros(4)
                for s in 1:4
                    Lab_guess[s]=Labor_solold[Zi,s]
                end
                opt = nlsolve(IM_FOC,Lab_guess,iterations=10000)
                if maximum(IM_FOC(opt.zero))>1e-6
                    println("bad zero, Zi=",Zi)
                #else
                #    println("good zero, Zi=",Zi)
                end

                Labor_sol[Zi,:]=opt.zero
                Vnew[Zi]=IM_Payoff(opt.zero,Z,Vold,Z_grid)
                λnew[Zi]=IM_λ(opt.zero,Z,λold,Z_grid)
                B_sol[Zi]=IM_B(opt.zero,Z,λold,Z_grid)
                nZ_sol[Zi,:]=IM_nZ(opt.zero,Z,λold,Z_grid)
            end
        end


        error=maximum(abs.(Vnew-Vold))
        println("iter ",iter,": error ",error)
        Vold,Vnew=Vnew,Vold
        λold,λnew=λnew,λold
    end

    policyfunctions=(λold,Labor_sol,B_sol,nZ_sol)

    return Vold,policyfunctions,failedpoints
end


###################################################################
# Forced Stationarity
function stat_solution(Labor,Z)

    error=zeros(4)
    cum=0
    for s in 1:4
        zs,gs=ExogStates(s)
        l=Labor[s]
        cons=zs*l-gs
        if cons<0
            error=1000*ones(4)
            return error
        end
        cum+=cons^(-α)*Pi[s]
    end

    for s in 1:4
        zs,gs=ExogStates(s)
        l=Labor[s]
        cons=zs*l-gs
        error[s]=(cons^(-α)/cum-β)*Z-(cons^(1-α)-l^(1+γ))
    end
    return error
end

function stat_payoff(Labor,Z)

    payoff=0
    for s in 1:4
        zs,gs=ExogStates(s)
        l=Labor[s]
        cons=zs*l-gs
        payoff+=cons^(1-α)/(1-α)-l^(1+γ)/(1+γ)*Pi[s]
    end
    payoff=payoff/(1-β)

    return payoff
end

function stat_consB(Labor,Z)

    cons=zeros(4)
    cum=0
    for s in 1:4
        zs,gs=ExogStates(s)
        l=Labor[s]
        cons[s]=zs*l-gs
        cum+=cons[s]^(-α)*Pi[s]
    end
    B=Z/cum

    return cons,B
end

#=

Labor_stat=zeros(n_Z,4)
Cons_stat=zeros(n_Z,4)
B_stat=zeros(n_Z)
V_stat=zeros(n_Z)

for Zi in 1:n_Z
    Z=New_Z_grid[Zi]
    f(x)=stat_solution(x,Z)
    Lab_guess=zeros(4)
    for s in 1:4
        if Zi==1
            lab_interp=LinearInterpolation(Sp_Results[:,1],Sp_Results[:,2+s],extrapolation_bc=Line())
            Lab_guess[s] = lab_interp(Z)
        else
            Lab_guess[s]=Labor_stat[Zi-1,s]
        end
    end
    opt=nlsolve(f,Lab_guess,iterations=5000)
    if maximum(f(opt.zero))>1e-6
        println("bad zero")
    end
    Labor_stat[Zi,:]=opt.zero
    Cons_stat[Zi,:],B_stat[Zi]=stat_consB(opt.zero,Z)
    V_stat[Zi]=stat_payoff(opt.zero,Z)
end

=#
