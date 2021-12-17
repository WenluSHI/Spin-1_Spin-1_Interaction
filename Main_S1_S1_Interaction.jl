### Beginning.
println('\n'); 
printstyled("======= Wenlu SHI Code ======="; color = :yellow); 
println('\n')

### License
include("section_Print_License.jl")

### Packages
using LinearAlgebra
using DelimitedFiles
include("section_Parameters.jl") # Parameters

### Define saving address
function func_save(x,y)	# value, name
	printstyled("Saving data ... "; color = :green);println(" ")
	address = join(["/Users/Lyint_Shi/Documents/JuliaLearning/Julia_s1_s1_interaction/",y,".tsv"])
	f=open(address,"w")
	writedlm(f,x,'\t')
	close(f)
	println("Success saving data!");
	println(address)
end

### Calculate
Points = 5000;
J = range(0,4,Points); #exchange interaction

S1p = sqrt(2)*diagm(1=>[1,1,0,1,1,0,1,1]);
S1m = sqrt(2)*diagm(-1=>[1,1,0,1,1,0,1,1]);
S1x = (S1p+S1m)/2;
S1y = (S1p-S1m)/2/im;
S1z = diagm(0=>[1,0,-1,1,0,-1,1,0,-1]);

s1p = diagm(3=>[1,1,1,1,1,1]);
s1m = diagm(-3=>[1,1,1,1,1,1]);
s1x = (s1p+s1m)/2;
s1y = (s1p-s1m)/2/im;
s1z = diagm(0=>[1,1,1,0,0,0,-1,-1,-1])./2;

# H = H_spin + J * H_interaction
H_spin = g1*mu_B*B*S1z + g2*mu_B*B*s1z + D1*(S1z*S1z) + D2*(s1z*s1z);
H_int = S1x*s1x+S1y*s1y+S1z*s1z;
H_spin = repeat(H_spin,1,1,Points);
H_int =  repeat(H_int,1,1,Points);
J = repeat(J,1,size(H_spin,1),size(H_spin,1));
J = permutedims(J,(3,2,1));
H = H_spin + J.* H_int;

#psi = zeros(size(H)); # eigen vectors
E = zeros(size(H,1),size(H,3)); # eigen energies

for j in 1:size(H,3)
	E[:,j] = eigvals(H[:,:,j]);
end

ee = zeros(size(E));
for j in 1:size(E,2)
	ee[:,j] = sort(E[:,j]);
end

Energy = ee';

### Save to file
data = [J[1,1,:] Energy];
func_save(data,"Energy_Level");
data = [J[1,1,:] Energy.-Energy[:,1]];
func_save(data,"Excitation_Energy");

### Plot
include("section_Plot.jl") 

### Clear
J, psi, data, ee, H_spin, H_int, H, Energy, E,= 0,0,0,0,0,0,0,0,0;
### Ending
println('\n');
printstyled("=========== Done ==========="; color = :yellow);
println('\n')



