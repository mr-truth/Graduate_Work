using QuadGK

na_array = [0.01, 0.001, 0.0001]
pump = 0.0625
compensation = 0.025
E_array = -4:0.12:8
t_array = [0, 10, 100, 1000]
PL = zeros(length(E_array))

for k=1:length(na_array)
na = na_array[k]
E_md = 2.515 * 2 * compensation * (na^(2/3))
println(na)
for j=1:length(t_array)
	for i=1:length(E_array)
		E = E_array[i]
		t = t_array[j]
		
		function integrate()
			func(energy_1) = (energy_1^-5) * (E_md/((4*(E_md/energy_1)^2) + (E - energy_1)^2)) * exp(-4 / energy_1 - 32 * π / 3 * na / energy_1 ^ 3 - t*exp(-4/energy_1)) 
			result = quadgk(func, 0, 1000, rtol=1e-12)[1]
			end
		
		function integrate_q(t)
			func(energy) = (exp(-t * exp(-4/energy)) - 1)/energy^4
			result = quadgk(func, 0, 1000, rtol=1e-12)[1]
			end
		
		a = integrate() * exp(integrate_q(t) * 4 * π *na)
		PL[i] = (64 * na * pump)/(32 * π * na) * a

		println(PL[i])
		end
	println("--------------")
	end
end