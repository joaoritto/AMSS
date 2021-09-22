# States: ll,lh,hl,hh = 1, 2, 3, 4

# we get our state variables as a function of s

function ExogStates(state)

    if state==1
        zs = z_l
        gs = g_l
    elseif state==2
        zs = z_l
        gs = g_h
    elseif state==3
        zs = z_h
        gs = g_l
    elseif state==4
        zs = z_h
        gs = g_h
    end

    return zs, gs

end

# We can express λ as function of l, and Z as a function of all the L's.
# Then we can get B as a function of l and Z.

function Lam_Lcm(labor,s)

    zs, gs = ExogStates(s)

    cons = zs*labor - gs

    if cons <= 0

        print("someone did something very bad here")

    end

    lam_1 =  cons^(-α)*zs - labor^γ

    lam_2 = (1-α)*cons^(-α)*zs - (1+γ)*labor^γ

    λ = -lam_1/lam_2

    return λ

end

function Z_Lcm(Labor)  # Labor is a 4-element vector for each state

    Z = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        if cons <= 0

            print("someone did something very bad here")

        end

        Z += ( cons^(1-α) - labor^(1+γ) )*Pi[s]/(1-β)

    end

    return Z

end

function Pay_Lcm(Labor)  # Labor is a 4-element vector for each state

    Payoff = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labor[s]

        cons = zs*labor - gs

        if cons <= 0

            print("someone did something very bad here")

        end

        Payoff += ( cons^(1-α)/(1-α) - labor^(1+γ)/(1+γ) )*Pi[s]/(1-β)

    end

    return Payoff

end

function Bond_LZcm(labor,Z,s)

    zs, gs = ExogStates(s)

    cons = zs*labor - gs

    if cons <= 0

        print("someone did something very bad here")

    end

    B = ( cons^(1- α) - labor^(1+γ) + β*Z ) / cons^(-α)

    return B

end

function Find_zero()

    # Finds approx values at Z = 0 using interpolation

    # We can solve for lambda given labor and s. We invert that to get our labor function

    # Size = 1000 # number of grid points

    Size=1000 # Increase this to have higher precision

    epi = .1 # shifter to prevent 0 consumption

    Results_Lam  = zeros(Size,4)

    Labors = zeros(Size,4)

    for state = 1:4

        zs, gs = ExogStates(state)

        for count = 1:Size

            facc = 1.7/Size  # Assuming a max labor of 2 so this scales our grid

            labor_in = facc*(count-1) + gs/zs + epi  # plausible labor levels after adjustment

            Results_Lam[count,state] = Lam_Lcm(labor_in,state)

            Labors[count,state] = labor_in

        end


    end


    function interp_LM(x,state)  # retuns labor in state as function of lambda

        X = reverse(Results_Lam[:,state]) # reversing order to make increasing

        Y = reverse(Labors[:,state])

        interp_cm = LinearInterpolation(X,Y)

        return interp_cm(x)

    end


    # Constructing the overall outcomes:

    # start with a counter over lambda values from which we get labors

    # use labors to get bond values

    # construct ex ante debt and Z values

    Labor = zeros(4)

    Range = 0:.0001:0.1

    ssize = size(Range)[1]

    Big_Res = zeros(ssize,10) # lambda + 4 labor values + 4 bond values + Z value


    for counter in 1:ssize  # start and stop choosen by looking at top and bottom above

        lam_count=Range[counter]
        Big_Res[counter,1] = lam_count

        for s = 1:4

            Labor[s] = interp_LM(lam_count,s)

            Big_Res[counter,s+1] = Labor[s]

        end

        Z = Z_Lcm(Labor)

        Big_Res[counter,10] = Z

        for a in 1:4
            Big_Res[counter,a+5] = Bond_LZcm(Labor[a],Z,a)
        end

    end



    # results with Z =0 and presumbably B as well

    find_0 = argmin(Big_Res[:,10].^2)  # trying to find the Z = 0 best approix.

    Case_0 = Big_Res[find_0,:]

    return Case_0,Big_Res

end

# Code for solving out given Guess and Z assuming constraint does not bind

# Here we define some functions to be used by the equation solver. This is to be used to
# make our solutions more exact, if necessary, and to check what we have.

function ERs_nb(Guess,Z)

    # This function generates the complete markets errors at guess,
    # given the state Z and the parameters para
    # Guess is a vector of labor levels and lambda

    Lambda = Guess[1]

    Labors = Guess[2:5]

    Errors = zeros(5)

    B = zeros(4)

    # The errors are with respect to the labor foc and implied Z

    Cum_Z = 0

    for s = 1:4

        zs, gs = ExogStates(s)

        labor = Labors[s]

        cons = zs*labor - gs

        Errors[s] = cons^(-α)*zs-labor^γ + Lambda*((1-α)*cons^(-α)*zs - (1+γ)labor^γ)

        B[s] = ( (cons^(1-α) - labor^(1+γ)) + β*Z ) / cons^(-α)

        Cum_Z += (cons^(1-α) - labor^(1+γ) + β*Z )*Pi[s]

    end

    Errors[5] = Cum_Z - Z

    return Errors

end

function Z_count(ZZ,debt_bound,Guess)

    # function constructs results counting up or down to
    # debt bound. ZZ is the grid which determines up or down


    ssize = size(ZZ)[1]

    sstuff = zeros(1,11)

    TempA = zeros(ssize,11)

    lastcount = 0

    for counter in 1:ssize+1

        zz = ZZ[counter]

        function F(x)  # function handle for fsolve

            return ERs_nb(x,zz)

        end

        results=nlsolve(F,Guess)

        answer = results.zero

        sstuff[1,1] = zz

        for count = 1:5
            sstuff[1,count+1] = answer[count]
        end

        bbonds = zeros(4)

        for s = 1:4

            bbonds[s] = Bond_LZcm(answer[s+1],zz,s)

            sstuff[1,6+s] = bbonds[s]

        end

        sstuff[11] = Pay_Lcm(answer[2:5])

        TempA[counter,:] = sstuff

        # working with positive and negative debt counds

        if debt_bound>0 && maximum(bbonds)>debt_bound


            lastcount = counter

            break

        elseif debt_bound<0 && minimum(bbonds)<debt_bound

            lastcount = counter

            break

        end

    end

    TempB = TempA[1:lastcount,:]

    return TempB

end


# Cases when Debt Bounds bind

# We need to interpolate μ=V'(\tilde Z)

function interp_MU(x)  # retuns μ  as a function of \tilde Z

    X = Sp_Results[:,1] # reversing order to make increasing

    Y = Sp_Results[:,2]

    interp_cm = LinearInterpolation(X,Y;extrapolation_bc=Line())

    return -interp_cm(x) # μ = - λ = V'(Z)

end


# Binds at the top of the Debt Range

function ERs_bbh(Guess,Z,bind_s)

    # This function generates the complete markets errors at guess,
    # given the state Z and the parameters para
    # Guess is a vector of lambda, labor levels and new Z.

    Lambda = Guess[1]

    Labors = Guess[2:5]

    newZ = Guess[6]

    Errors = zeros(6)

    B = zeros(4)

    # The errors are with respect to the labor foc and implied Z

    Cum_Z = 0

    for s = 1:4  # exceeds bound in ExogStates 4

        zs, gs = ExogStates(s)

        labor = Labors[s]

        cons = zs*labor - gs

        if s == bind_s

            Mu = interp_MU(newZ)  # V'(Z) = mu

            Cum_Z += (cons^(1-α) - labor^(1+γ) + β*newZ )*Pi[s]

            fac1 = (Mu+Lambda)*α*cons^(-α-1)*zs*M_max

            fac2 = Mu*((1-α)*cons^(-α)*zs - (1+γ)labor^γ)

            Errors[s] = cons^(-α)*zs-labor^γ - fac1 - fac2

            B[s] = M_max

            Errors[6] = cons^(-α)*M_max - cons^(1-α) + labor^(1+γ) - β*newZ

        else

            Errors[s] = cons^(-α)*zs-labor^γ + Lambda*((1-α)*cons^(-α)*zs - (1+γ)labor^γ)

            B[s] = ( (cons^(1-α) - labor^(1+γ)) + β*Z ) / cons^(-α)

            Cum_Z += (cons^(1-α) - labor^(1+γ) + β*Z )*Pi[s]

        end

    end

    Errors[5] = Cum_Z - Z

    return Errors

end


# working on the bottom end of the bounds

function ERs_bbl(Guess,Z,bind_s)

    # This function generates the complete markets errors at guess,
    # given the state Z and the parameters para
    # Guess is a vector of lambda, labor levels and new Z.

    Lambda = Guess[1]

    Labors = Guess[2:5]

    newZ = Guess[6]

    Errors = zeros(6)

    B = zeros(4)

    # The errors are with respect to the labor foc and implied Z

    Cum_Z = 0

    for s = 1:4  # exceeds bound in state 4

        zs, gs = ExogStates(s)

        labor = Labors[s]

        cons = zs*labor - gs

        if s == bind_s

            Mu = interp_MU(newZ)  # V'(Z) = mu

            Cum_Z += (cons^(1-α) - labor^(1+γ) + β*newZ )*Pi[s]

            fac1 = (Mu+Lambda)*α*cons^(-α-1)*zs*M_min

            fac2 = Mu*((1-α)*cons^(-α)*zs - (1+γ)labor^γ)

            Errors[s] = cons^(-α)*zs-labor^γ - fac1 - fac2

            B[s] = M_min

            Errors[6] = cons^(-α)*M_min - cons^(1-α) + labor^(1+γ) - β*newZ

        else

            Errors[s] = cons^(-α)*zs-labor^γ + Lambda*((1-α)*cons^(-α)*zs - (1+γ)labor^γ)

            B[s] = ( (cons^(1-α) - labor^(1+γ)) + β*Z ) / cons^(-α)

            Cum_Z += (cons^(1-α) - labor^(1+γ) + β*Z )*Pi[s]

        end

    end

    Errors[5] = Cum_Z - Z

    return Errors

end
