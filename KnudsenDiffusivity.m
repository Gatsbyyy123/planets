syms D_K D_Ki D_Kj %Knudsen diffusivity, of species i, and of spcies j (m^2/s)
syms R_gT
function D_k=Knudsen_diffusivity(r)
D_ki=4/3*((8*R_gT)/(pi*MW_i))^(1/2)*((3*pi*r)/(4*(pi+8)*(1-epsilon))) %Knudsen diffusivity of species j
D_kj=4/3*((8*R_gT)/(pi*MW_j))^(1/2)*((3*pi*r)/(4*(pi+8)*(1-epsilon))) %Knudsen diffusivity of species j

D_k=1/D_ki+1/D_kj
end