begin #space discretization
  L = 1
  Nx = 100
  xl = LinRange(0,L,Nx+2)
  xlu = xl[2:end-1]
  Deltax = L/(Nx+1)
end
#comentário test
