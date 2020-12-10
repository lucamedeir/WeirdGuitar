begin #space discretization
  L = 1
  Nx = 100
  xl = LinRange(0,L,Nx+2)
  xlu = xl[2:end-1]
  Deltax = L/(Nx+1)
end


begin #time discretization
    T=10
    Nt=100
    tl=LinRange(0,T,Nt)
    Δt=T/(Nt-1)
end;

c(x) =17; #velocity profile

begin #discretization matrix
    usingLinearAlgebra
    M11=zeros(Nx,Nx)
    M12=IM21=Tridiagonal((c.(xlu[2:end]).^2)./Δx^2,(c.(xlu).^2).*(-2)/Δx^2,(c.(xlu[1:end-1]).^2)./Δx^2)
    M22=zeros(Nx,Nx
    )M= [M11M12;M21M22]
end;

begin #implicit Euler
    methodMM=I+M.*(Δt)
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

gif(anim,"Ruela's_guitar.gif",fps=10)
