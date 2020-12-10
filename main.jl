begin #space discretization
  L = 1
  Nx = 100
  xl = LinRange(0,L,Nx+2)
  xlu = xl[2:end-1]
  Δx = L/(Nx+1)
end

begin #time discretization
    T=1
    Nt=100
    tl=LinRange(0,T,Nt)
    Δt=T/(Nt-1)
end;

c(x) =17; #velocity profile

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

begin #implicit Euler method
    MM=I+M.*(Δt)
    iMM=inv(Array(MM))
end;

begin #initial conditions
    phi0=sin.((pi/L).*xlu)
    pi0=zeros(Nx)
    sol= [[phi0; pi0]]
end;

begin #applies implicit Euler method
    for i in 1:Nt
        push!(sol,iMM*sol[i][:])
    end
end;

begin #animation
    using Plots
    anim=Plots.Animation()
    for i in 1:Nt
        plot(xl,[0;sol[i][1:Nx];0],legend=false,ylims=(-1,1))
        Plots.frame(anim)
    end
end

gif(anim,"Ruela's_guitar.gif",fps=30)
