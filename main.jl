#===============================================================================
NUMERICAL SETUP
===============================================================================#
begin # Space discretization
  L = 1 # Length of unstretched string (in meters)
  Nx = 1001 # Nx must be an odd number if you want to use Simpson's rule to calculate the total "energy" later (as it is currently being done)
  xl = LinRange(0,L,Nx+2) # Space grid including boundary points
  xlu = xl[2:end-1] # Space grid without boundary points
  Δx = L/(Nx+1) # Space step
end;

begin # Time discretization
    T = 10 # Maximum elapsed time (in seconds)
    Nt = 1000
    tl = LinRange(0,T,Nt) # Time grid
    Δt = T/(Nt-1) # Time step
end;

begin # Speed of sound function
    function c(x) # Can be changed at will
        1/10
    end
end;
#===============================================================================
SOLVING THE PROBLEM
===============================================================================#
begin # Building the discretization matrix M (see LaTeX document)
    using LinearAlgebra
    M11 = zeros(Nx,Nx)
    M12 = I
    beta1 = (c.(xlu[2:end]).^2)./Δx^2
    alfa = (c.(xlu).^2).*(-2)/Δx^2
    beta2 = (c.(xlu[1:end-1]).^2)./Δx^2
    M21 = Tridiagonal(beta1,alfa,beta2)
    M22 = zeros(Nx,Nx)
    M = [M11 M12;M21 M22]
end;

begin # Crank-Nicolson method
    MMl = I-M.*(Δt/2) # Left-hand side matrix
    MMr = I+M.*(Δt/2) # Right-hand side matrix
    iMM = inv(Array(MMl))*MMr
end;

begin # Initial conditions (can be chosen at will)
    function pluck(x) # Initial triangular-shaped string
        px = 2*L/3 # Position of the vertex
        py = 1 # Height of the vertex
        if x < px
            py/px*x
        else
            py/(L-px)*(-x+L)
        end
    end
    function gauss(x) # Initial Gaussian-shaped string
        μ = L/2 # Mean
        σ = L/100 # Standard deviation
        exp(-(x-μ)^2/(2*σ)^2)
    end
    function dgauss(x) # Time derivate of gauss(x)
        μ = L/2 # Mean
        σ = L/100 # Standard deviation
        exp(-(x-μ)^2/(2*σ)^2)*(-2)/(2*σ^2)*(x-μ)*(c(x))
    end
    phi0 = sin.((pi/L).*xlu)#gauss.(xlu)#=pluck.(xlu)=##=# # Initial position Φ0
    pi0 = #=dgauss.(xlu)=#zeros(Nx) # Initial velocity Π0
    sol = [[phi0; pi0]] # Initial Y vector (see LaTeX document)
end;

begin # Time evolution
    for i in 1:(Nt-1)
        push!(sol,iMM*sol[i][:])
    end
end;
#===============================================================================
CALCULATING THE "ENERGY" FUNCTION (SANITY CHECK)
===============================================================================#
begin # Separating the position and velocity parts from the Y vector
    PHI = [] # Position Φ
    PI = [] # Velocity Π
    for i in 1:Nt
        push!(PHI,[0;sol[i][1:Nx];0])
        push!(PI,[0;sol[i][Nx+1:end];0])
    end
end

begin # Calculating the gradient ∇Φ for all time steps (excluding the boundaries)
    dPHI = [] # Gradient ∇Φ
    for i in 1:Nt
        aux = (PHI[i][3:end].-PHI[i][1:end-2])./(2*Δx) # Second-order formula
        push!(dPHI,aux)
    end
end;

begin # Calculating "energy" density function
    energy = [] # "Energy" density
    for i in 1:Nt
        push!(energy,(1)./((c.(xlu)).^2).*((PI[i][2:end-1]).^2).+((dPHI[i][:]).^2)) # See LaTeX document
    end
end

begin # Simpson's rule for total "energy"
    ENERGY = [] # Integral of "energy" density for all time steps
    Simpson = [1;zeros(Nx-2);1] # Vector of weights from Simpson's formula
    for i in 2:Nx-1
        if i % 2 == 0
            Simpson[i] = 4
        else
            Simpson[i] = 2
        end
    end
    Simpson = Simpson.*(Δx/3)
    for i in 1:Nt
        push!(ENERGY,dot(Simpson,energy[i]))
    end
end
#===============================================================================
PLOTTING THE RESULTS
===============================================================================#
begin # Plotting total "energy" X time (sanity check)
    using Plots
    plot(tl,ENERGY,label="Total Energy")
end

begin # Plotting string and "energy" density together
    anim = Plots.Animation()
    for i in 1:Nt
        plot(xlu,PHI[i][2:end-1],label="Guitar string",ylims=(-1,1)) # String
        plot!(xlu,(energy[i][:].-(ENERGY[1]/L))./maximum(energy[1][:]),label="Energy density") # "Energy" density (rescaled)
        Plots.frame(anim)
    end
end

gif(anim,"Ruela's_guitar.gif",fps#==Nt/T=#=30)
