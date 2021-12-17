printstyled("Making plots ... "; color = :green);println(" ")
##################################################################
using PyPlot

figure(1)
plot(J[1,1,:],Energy)
legend(["a","ab","abc"]);
xlabel("Exchange interaction J (meV)");
ylabel("Energy (meV)")

figure(2)
plot(J[1,1,:],Energy.-Energy[:,1])
xlabel("Exchange interaction J (meV)");
ylabel("Excitation Energy (meV)")
##################################################################
println("Success making plots!");