begin #space discretization
  L = 1
  Nx = 1000
  xl = LinRange(0,L,Nx+2)
  xlu = xl[2:end-1]
  Δx = L/(Nx+1)
end

begin #time discretization
    T=1
    Nt=500
    tl=LinRange(0,T,Nt)
    Δt=T/(Nt-1)
end;

c(x) =1; #velocity profile

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
    phi0=sin.(( 50 * pi/L).*xlu) + sin.(( 123 * pi/L).*xlu)
    pi0=zeros(Nx)
    sol= [[phi0; pi0]]
end;

begin #applies implicit Euler method
    for i in 1:Nt
        push!(sol,iMM*sol[i][:])
    end
end;

begin
    using FFTW
    fft_sol = []
    for i in 1:Nt
        push!(fft_sol, fft([0;sol[i][1:Nx];0]) |> fftshift)
    end

end

begin #animation
    using Plots
    anim=Plots.Animation()

    Xs = 1 / (1.1 * (Nx+2)/2)
    freqs = fftfreq(length(xl), 1.0/Xs) |> fftshift

    for i in 1:Nt
        plot(freqs,abs.(fft_sol[i]),legend=false,xlim=(0,200),ylim=(0,500))
        Plots.frame(anim)
    end
end

gif(anim,"Ruela's_guitar.gif",fps=30)
