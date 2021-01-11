begin #space discretization
  L = 1
  Nx = 1000
  xl = LinRange(0,L,Nx+2)
  xlu = xl[2:end-1]
  Δx = L/(Nx+1)
end;

begin #time discretization
    T=10
    Nt=100
    tl=LinRange(0,T,Nt)
    Δt=T/(Nt-1)
end;

begin
    function c(x) #speed of sound
        1/10
    end
end;

begin #discretization matrix
    using LinearAlgebra
    M11=zeros(Nx,Nx)
    M12=I
    beta1 = (c.(xlu[2:end]).^2)./Δx^2
    alfa = (c.(xlu).^2).*(-2)/Δx^2
    beta2 = (c.(xlu[1:end-1]).^2)./Δx^2
    M21=Tridiagonal(beta1,alfa,beta2)
    M22=zeros(Nx,Nx)
    M= [M11 M12;M21 M22]
end;

begin #Crank-Nicolson method
    MMl=I-M.*(Δt/2) #left-hand side matrix
    MMr=I+M.*(Δt/2) #right-hand side matrix
    iMM=inv(Array(MMl))*MMr
end;

begin #initial conditions
    function guitar(x)
        px = 2*L/3
        py = 1
        if x < px
            py/px*x
        else
            py/(L-px)*(-x+L)
        end
    end
    function gauss(x)
        μ = L/2
        σ = L/100
        exp(-(x-μ)^2/(2*σ)^2)
    end
    function dgauss(x)
        μ = L/2
        σ = L/100
        exp(-(x-μ)^2/(2*σ)^2)*(-2)/(2*σ^2)*(x-μ)*(c(x))
    end
    phi0=gauss.(xlu)#=guitar.(xlu)=##=sin.((pi/L).*xlu)=#
    pi0=#=dgauss.(xlu)=#zeros(Nx)
    sol= [[phi0; pi0]]
end;

begin #applies implicit Euler method
    for i in 1:(Nt-1)
        push!(sol,iMM*sol[i][:])
    end
end;

begin #separate the string position and its velocity and include the boundary values
    PHI = []
    PI = []
    for i in 1:Nt
        push!(PHI,[0;sol[i][1:Nx];0])
        push!(PI,[0;sol[i][Nx+1:end];0])
    end
end

begin #calculates the gradient for all time steps (second order inside the domain)
    dPHI = []
    for i in 1:Nt
        aux = (PHI[i][3:end].+PHI[i][1:end-2])./(2*Δx)
        push!(dPHI,aux)
    end
end;

begin #calculates "energy" function
    energy = []
    for i in 1:Nt
        push!(energy,((1)./((c.(xlu)).^2).*((PI[i][2:end-1]).^2).+((dPHI[i][:]).^2)))
    end
end

#=begin
    using Plots
    Eanim=Plots.Animation()
    for i in 1:Nt
        plot(xlu,energy[i][:]./100000,legend=false,ylims=(-0.5,10))
        Plots.frame(Eanim)
    end
end=#

begin #animation
    using Plots
    anim=Plots.Animation()
    for i in 1:Nt
        plot(xlu,PHI[i][2:end-1],label="Guitar string",ylims=(-1,1))
        plot!(xlu,energy[i][:]./1000000,label="Energy density")
        Plots.frame(anim)
    end
end

#gif(Eanim,"Ruela's_guitar_energy.gif",fps=Nt/T#=30=#)
gif(anim,"Ruela's_guitar.gif",fps=Nt/T#=30=#)
