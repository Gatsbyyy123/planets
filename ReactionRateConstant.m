syms A Ea %frequency factor (/s), apparant activation energy (J/mol)
k=zeros(1,10)

for i=1:10
    k(i)=reaction_rate_constant(T)
end
ad=zeros(4)
ad(3,4)=3
function kt=reaction_rate_constant(T)
A=40
Ea=72
T=873.15
R_gT=8.3145
kt=A*exp(-Ea/(R_gT*T))
end
