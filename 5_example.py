from chapters.ffivl import critical_pt_viscocity
class T:
	def test_2_1(s):
		# Calculate the viscocity of nitrogen of 20°C and amtospheric pressure from 2.8
		P_c = 33.5
		T_c = 126.2 #Kelvin
		M = 28.0 #Nitrogen
		ans = critical_pt_viscocity(M=M, P_c=P_c, T_c=T_c)
		print (ans)

if __name__=="__main__":
	T().test_2_1()
