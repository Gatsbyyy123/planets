syms T P v_i v_j MW_i MW_j %temperature (K), system pressure (Psi), volume fractions of species i (vol.%), volume fractions of species j (vol.%), molecualr mass of i (kg/mol), and j (kg/mol)
syms D_i D_j D_ij %diffusivity of single component and binary diffusivity
syms n_s %number of species
syms x_i x_j %mole fractions of species i and j (mol.%)
syms V_i V_j %volume of species i and j (m^3)
syms n_i n_j %mole number of species i and j (mol)
syms rho_i rho_j %density of gas species (kg/m^3)
function D_ij = binary_diffusivity(T,MW_i,MW_j,P,n_i,n_j)
D_ij = 1.01325*10^2*(T^1.75*(1/MW_i +1/MW_j )^(1/2))/(P*(v_i^(1/3)+v_j^(1/3))^2)
D_i = 1/(1/(1-x_i)*(x_i/D_ij))
D_j = 1/(1/(1-x_j)*(x_i/D_ij))
%assume that the gas species are ideaAl gases
V_i=(n_i*R_gT*T)/P
V_j=(n_j*R_gT*T)/P
v_i =V_i/(V_i+V_j)
v_j =V_j/(V_i+V_j)
%n_i=rho_i*vi*(V_i+V_j)/WM_i
%n_j=rho_j*vj*(V_i+V_j)/WM_j
x_i=n_i/(n_i+n_j)
x_j=n_j/(n_i+n_j)
end