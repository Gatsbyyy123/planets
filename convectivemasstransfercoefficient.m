syms N_sh N_Re N_Sc % Sherwood number, Reynolds number, and Schmidt number
syms Velocity L % velocity of gas flow (m/s), and characteristic length of solid particle (m)
syms kvis_i % kinetic viscosity of species i (m^2/s)
syms K_m %convective mass transfer coefficient (m/s)
function K_m = convective_mass_transfer_coefficient(Velocity,kvis_i)
N_sh=2.0+0.6*N_Re^(1/2)*N_Sc^(1/3)
N_Re=Velocity*L/(kvis_i)
N_Sc=(kvis_i)/D_ij
K_m=N_sh*D_ij/L
end