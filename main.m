syms T P v_i v_j MW_i MW_j %temperature (K), system pressure (Psi), volume fractions of species i (vol.%), volume fractions of species j (vol.%), molecualr mass of i (kg/mol), and j (kg/mol)
syms D_i D_j D_ij %diffusivity of single component and binary diffusivity
syms n_s %number of species
syms x_i x_j %mole fractions of species i and j (mol.%)
syms V_i V_j %volume of species i and j (m^3)
syms n_i n_j %mole number of species i and j (mol)
syms rho_i rho_j %density of gas species (kg/m^3)

particle_diameter=5, particle_radius=particle_diameter/2 %unit [mm]
simulation_time=120 %unit [s]
porosity=0.2, tortuosity=2, kvis_i=0.88*10^(-5)
bulk_pressure=101325, fluid_flow=200, temeprature=1073.13 %unit [Pa], unit [mL/min], %unit [K]
MW_FeO=72, MW_H2=2, MW_Fe=56, MW_H2O=18 %unit [g/mol]
pho_FeO=5.88*10^6, pho_H2=80.77 %unit [g/m^3]
alpha=1/3 %characteristic parameter
stoichio_H2=1, stoichio_H2O=1, stoichio_FeO=1, stoichio_Fe=1 %stoichiometric number
volume=4/3*pi*(particle_radius/1000)^3 %unit [m^3]
surface_area=4*pi*(particle_radius/1000)^2 %unit [m^2]
mass_FeO=(1-porosity)*volume*pho_FeO %unit [g]
initial_fraction_H2=0.999, initial_fraction_H2O=1-initial_fraction_H2
initial_n_H2=0.01, initial_n_FeO=mass_FeO/MW_FeO, initial_n_H2O=(1-initial_fraction_H2)*initial_n_H2 %unit [mol]
R_gT=8.3145 % ideal gas constant (J*K^-1*mol^-1)
C_H2= bulk_pressure/(R_gT*temeprature) %assumption: hydrogen is ideal gas
C_FeO=initial_n_FeO/volume %unit [mol/m^3]

time_counter=0
%transport phenomena properties
if time_counter==0
    n_H2=initial_n_H2O
    n_FeO=initial_n_FeO
    n_H2O=initial_n_H2O
end 

D_H2_H2O=binary_diffusivity(temeprature,MW_H2,MW_H2O,bulk_pressure,n_H2,n_H2O)
convective_mass_transfer_co=convective_mass_transfer_coefficient(N_Re,N_sc)


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