format shortE
particle_diameter=2; particle_radius=particle_diameter/2; %unit [mm]
total_simulation_time=30; %unit [s]
porosity=0.2; tortuosity=2; kinematic_viscosity_H2=1.84*10^(-5); L=particle_diameter/1000; %unit [Pa*s=N*s/m^2], unit[m]
bulk_pressure=101325; fluid_flow=200; temperature=1073.13; %unit [Pa], unit [mL/min], %unit [K]
cross_sectional_area=pi*(particle_radius/1000)^2;
bluk_flow_velocity=fluid_flow/60/10^6/cross_sectional_area;
MW_FeO=72; MW_H2=2; MW_Fe=56; MW_H2O=18; %unit [g/mol]
pho_FeO=5.88*10^6; pho_H2=80.77; %unit [g/m^3]
alpha=1/3; %characteristic parameter
A=40; Ea=72000; % reaciton kinetics unit [1/s] unit [J/mol]
stoichio_H2=1; stoichio_H2O=1; stoichio_FeO=1; stoichio_Fe=1; %stoichiometric number
volume=4/3*pi*(particle_radius/1000)^3; %unit [m^3]
surface_area=4*pi*(particle_radius/1000)^2; %unit [m^2]
mass_FeO=(1-porosity)*volume*pho_FeO; %unit [g]
initial_fraction_H2=0.999; initial_fraction_H2O=1-initial_fraction_H2;
initial_n_H2=0.01; initial_n_FeO=mass_FeO/MW_FeO; initial_n_H2O=(1-initial_fraction_H2)*initial_n_H2; %unit [mol]
R_gT=8.3145; % ideal gas constant (J*K^-1*mol^-1)
C_H2= bulk_pressure/(R_gT*temperature); %assumption: hydrogen is ideal gas
C_FeO=initial_n_FeO/volume; %unit [mol/m^3]


%%%transport phenomena properties
D_H2_H2O=binary_diffusivity(temperature,MW_H2,MW_H2O,bulk_pressure,initial_n_H2,initial_n_H2O,R_gT)
K_m=convective_mass_transfer_coefficient(bluk_flow_velocity,kinematic_viscosity_H2,L,D_H2_H2O)
D_k=Knudsen_diffusivity(MW_H2,MW_H2O,particle_radius,porosity,R_gT)
D_overall=overall_diffusivity(temperature,MW_H2,MW_H2O,bulk_pressure,D_k,initial_n_H2,initial_n_H2O,R_gT) %ALERT: the values of diffusivity are obstrusively high

%%%settings of the time steps and mesh size, as well as the acquisition
%rate.
number_elements=10;
delta_r=particle_radius/number_elements %meaning that on the radial direction, the interval is divided into 50 segments.
total_time_steps=ceil((2*D_overall*total_simulation_time)/(porosity*delta_r^2)) %stability criterion
delta_t=total_simulation_time/total_time_steps
element_location=zeros(number_elements+1) %???????? +1 or not?
for f=1:1:number_elements
    element_location(1)=delta_r;
    element_location(f)=element_location(f)+delta_r;
end 


%%%the following part is to creat matrice for data storage, and the
%calculated values at every 1000 time steps will be returned, because the
%total time steps are way too large. This can help to limit the size of
%matrices.
mass_transfer_rate_array=zeros(ceil(total_time_steps/100),1);
mole_number_FeO_array=zeros(ceil(total_time_steps/100),1);
consumed_mole_FeO_array=zeros(ceil(total_time_steps/100),1);
mass_loss_percentage_array=zeros(ceil(total_time_steps/100),1);
average_conversion_degree_array=zeros(ceil(total_time_steps/100),1);
consumption_rate_FeO_array=zeros(ceil(total_time_steps/100),1);
average_reaction_rate_array=zeros(ceil(total_time_steps/100),1);
plotting_exporting_time_array=zeros(ceil(total_time_steps/100),1);


%%%IMPORTANT every column of the following two matrice represents a mesh
%element for seperate data storage at each time step.
FeO_concentration_mesh_element_array=ones(ceil(total_time_steps/100),number_elements);
H2_concentration_mesh_element_array=ones(ceil(total_time_steps/100),number_elements);

%%%calculations of the initial values

initial_mass_FeO=mass_FeO;
initial_concentration_H2=C_H2;
initial_concentration_FeO=C_FeO;
max_mass_loss_percentage=(initial_mass_FeO-initial_n_H2/stoichio_FeO*stoichio_Fe*MW_Fe)/initial_mass_FeO;
concentration_H2=zeros(number_elements,1); %hydrogen was afterwards diffused into the particle, and therefore its initial concentration in the particle is zero, while the FeO is opposite.
concentration_FeO=ones(number_elements,1);
concentration_Fe=ones(number_elements,1);
concentration_FeO=concentration_FeO*initial_concentration_FeO;
concentration_Fe=initial_concentration_FeO*concentration_Fe-concentration_FeO;

reaction_rate_array=zeros(ceil(total_time_steps/100),1); %??????????(revised) every time step should have a reaction rate.
%consumption_rate_FeO_array=zeros(ceil(total_time_steps/100),1); % why to initialise these values at the beginning of each time step? (revised)

time_step_counter=0;
element_counter=0;
simulation_time=0;
%%%time looping starts from now. The following calculation represent the
%data acquirement of all mesh elements in a single time step.
for i=1:1:total_time_steps
 
    if time_step_counter==0 %initial values
    n_H2=initial_n_H2;
    n_FeO=initial_n_FeO;
    n_H2O=initial_n_H2O;
    n_Fe=0;
    end 

    if time_step_counter==total_time_steps
        disp=(['the final step is ' num2str(time_step_counter), ' and the simulation is ' num2str(time_step_counter/total_time_steps*100) '% completed'])
    end 

time_step_counter=time_step_counter+1
simulation_time=simulation_time+delta_t
element_counter=0;

% to calculate the reaction rate for the current time step, from the values
% 
    

for j=(number_elements):-1:1 
    %this is to solve the values of each element, including reaction rate, concentrations at the surface, 
    % at each mesh element (by gas species mass balance + central difference method)
    element_counter=element_counter+1

    %conversion_degree=1-(concentration_FeO(j)/initial_concentration_FeO)
    %f_x=(1-conversion_degree)^(1/3) %used as the shrinking core model
    
    k=reaction_rate_constant(A,Ea,temperature,R_gT);
    
    if j>=(number_elements-1)   %surface and subsurface
    concentration_H2(number_elements)=initial_concentration_H2;
    f_x=(concentration_FeO(number_elements)/initial_concentration_FeO)^(1/3);
    reaction_rate_array(number_elements)=stoichio_H2*k*concentration_H2(number_elements)*f_x ; %reaction rate is never zeron on surface
    
    %B.C. on the surface, includes the 2nd-order backward difference of K_m&
    %2nd-order central difference method of D_overall.
    f_x=(concentration_FeO(number_elements-1)/initial_concentration_FeO)^(1/3);
    M=(element_location(number_elements)^2+element_location(number_elements-1)^2)/element_location(number_elements)^2;
    N=stoichio_FeO/stoichio_H2*k*f_x/(K_m/delta_r-D_overall/delta_r^2);
    P=concentration_H2(number_elements)+element_location(number_elements-1)^2/element_location(number_elements)^2*concentration_H2(number_elements-2);
    concentration_H2(number_elements-1)=P/(M-N)
    
    reaction_rate_array(number_elements-1)=stoichio_H2*k*concentration_H2(number_elements-1)*f_x
    
    %the hydroge concentration at the subsurface element is still zero and waiting to be updated by Gas Phase Mass Balance.
    %boundary conditions of mass transfer and diffusion, to obtain the
    %concentration of specific hydrogen at the surface element.
    a_0=4/3*pi*(element_location(j)/1000)^3/(4*pi*(element_location(j)/1000)^2); %the geometry characteristic parameter

    if (stoichio_FeO/stoichio_H2*(1-porosity)*reaction_rate_array(j)*delta_t*1/a_0)>= initial_concentration_FeO % FUSE -- if the chemical consumption on surface is too high
        concentration_FeO(j)=0      
    else 
        concentration_FeO(j)=concentration_FeO(j)-stoichio_FeO/stoichio_H2*(1-porosity)*reaction_rate_array(j)*delta_t*1/a_0
    end 

    end

    if (j>=2)&&(j<=number_elements-2)  %middle
        a_1=4/3*pi*(element_location(j)/1000)^3/(4*pi*(element_location(j)/1000)^2)
                
        B=delta_t/porosity*D_overall/(element_location(j+1)^2)
        A=1+B*element_location(j)^2/(2*delta_r^2)
        concentration_H2(j)=((1+delta_t/porosity*porosity/a_1*k)*concentration_H2(j+1)-B*((concentration_H2(j+2)-concentration_H2(j+1))*element_location(j+2)^2)/(2*delta_r^2)+B*element_location(j)^2*concentration_H2(j+1)/(2*delta_r^2))/A
        %mass balance of gas species inside the particle, and therefore the
        %convective mass transfer is excluded. remember to consider the
        %porosity. Additionally, the calculation should start from the 2nd
        %element because of the central B.C., so it starts from n+1.
        reaction_rate_array(j)=stoichio_H2*k*concentration_H2(j)*f_x %update the reaction rate at the new element
        
    if (stoichio_FeO/stoichio_H2*(1-porosity)*reaction_rate_array(j)*1/a_1*delta_t)>=concentration_FeO(j) %FUSE
        concentration_FeO(j)=0
    else 
        concentration_FeO(j)=concentration_FeO(j)-(1-porosity)*reaction_rate_array(j)*stoichio_FeO/stoichio_H2*1/a_1*delta_t
    % consumption of FeO due to chemical reaction.
    end       
    end

    if j==1 % Subcentre and centre
    a_2=4/3*pi*(element_location(j)/1000)^3/(4*pi*(element_location(j)/1000)^2);
    concentration_H2(1)=concentration_H2(2) % B.C. the concentration gradient at the centre is zero.
    reaction_rate_array(j)=stoichio_H2*k*concentration_H2(j)*f_x
    if (stoichio_FeO/stoichio_H2*(1-porosity)*reaction_rate_array(j)*1/a_2*delta_t)>=concentration_FeO(j) %FUSE
        concentration_FeO(j)=0
    else 
        concentration_FeO(j)=concentration_FeO(j)-(1-porosity)*reaction_rate_array(j)*stoichio_FeO/stoichio_H2*1/a_2*delta_t
    end
    end
    
    H2_concentration_mesh_element_array(j,i)=concentration_H2(j) 
    FeO_concentration_mesh_element_array(j,i)=concentration_FeO(j)
    %This means store the value of j^th mesh element at the i^th time step
    %and j^th mesh element.Therefore, as time marching, the values would be
    %stored at the next row, and the new time step would run to fill every
    %mesh element with new a result.
end

%Data organising
mole_number_FeO=0; % mole nuumber of FeO should be reset at every time step to pass down the data for storage and output.
for m=(number_elements-1):-1:2
        mole_number_FeO=(4/3*pi*(element_location(m)/1000)^3-4/3*pi*(element_location(m-1)/1000)^3)*concentration_FeO(m)+mole_number_FeO
        mole_consumption_rate=(4/3*pi*(element_location(m)/1000)^3-4/3*pi*(element_location(m-1)/1000)^3)*reaction_rate_array(m)*(1-porosity)/stoichio_H2*stoichio_FeO
    
    average_mole_number=mole_number_FeO/(volume)
    consumed_mole_FeO=initial_n_FeO-mole_number_FeO %scalar, which would be used further for data storage in matrice.
    average_conversion_degree=consumed_mole_FeO/initial_n_FeO
    mass_loss=consumed_mole_FeO*MW_FeO  %unit [g]
    mass_loss_percentage=mass_loss/initial_mass_FeO
    average_reaction_rate=mole_consumption_rate/volume 
end 
   
% transfer the results to storage matrice at the current time step.
if rem(i,100)==0||i==total_time_steps
    acquisition_index=int16(i/100);
    mass_transfer_rate_array(acquisition_index)=surface_area*K_m*(initial_n_H2-concentration_H2(number_elements-1));
    mole_number_FeO_array(acquisition_index)=mole_number_FeO;
    consumed_mole_FeO_array(acquisition_index)=consumed_mole_FeO;
    mass_loss_percentage_array(acquisition_index)=mass_loss_percentage;
    average_conversion_degree_array(acquisition_index)=average_conversion_degree;
    consumption_rate_FeO_array(acquisition_index)=mole_consumption_rate;
    average_reaction_rate_array(acquisition_index)=average_reaction_rate;
    plotting_exporting_time_array(acquisition_index)=i;
    if rem(acquisition_index,10)==0  %every 1000 time steps, output all the data in workspace.
       filename = 'output_data';
       headers={'mass_transfer_rate_array [mol/s]','mole_number_FeO_array [mol]','consumed_mole_FeO_array [mol]','mass_loss_percentage_array [wt.%]','average_conversion_degree_array [%]','consumption_rate_FeO_array [mol/s]','average_reaction_rate_array [mol/(m^3*s)]'};
       compile_data=[mass_transfer_rate_array,mole_number_FeO_array,consumed_mole_FeO_array,mass_loss_percentage_array,average_conversion_degree_array,consumption_rate_FeO_array,average_reaction_rate_array];
       combination=strcat(filename,num2str(round(acquisition_index/10)),'.xlsx');
       mkdir results
       writecell(headers,combination,'Sheet',1,'Range','A1');
       writematrix(compile_data,combination,'Sheet',1,'Range','A2');
       movefile output_data*.xlsx results\

       save(fullfile('results/', ['workspace_' num2str(round(acquisition_index/10))]));
    end 
end

end 
