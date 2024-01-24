using QuadGK

na = 0.00235
pump = 0.99
compensation = 0.17647058823529413
E_md = 2.515 * 2 * compensation * (na^(2/3))
E_array = -4:0.12:8
t_array = [5, 15, 50, 500, 5000, 50000, 500000]
PL = zeros(length(E_array))

for j=1:length(t_array)
	for i=1:length(E_array)
		E = E_array[i]
		t = t_array[j]
		
		function integrate_q()
			func(energy) = (exp(-t * exp(-4/energy)) - 1)/energy^4
			result = quadgk(func, 0, 1000, rtol=1e-12)[1]
			end
		
		function integrate()
			func(energy_1) = (energy_1^-5) * ((E_md*(1-integrate_q()))/((4*((E_md*(1-integrate_q()))/energy_1)^2) + (E - energy_1)^2)) * exp(-4 / energy_1 - 32 * π / 3 * na / energy_1 ^ 3 - t*exp(-4/energy_1)) 
			result = quadgk(func, BigFloat(0), 1000, rtol=1e-12)[1]
			end
		
		a = integrate() * exp(integrate_q() * 4 * π *na)
		PL[i] = (64 * na * pump)/(32 * π * na) * a

		println(PL[i])
		end
	println("--------------")
	end
