d = {10:21, 11:5, 9:6, 8:8, 7:18, 6:11, 5:16, 4:11, 3:17}

def ter(a,i):
	return f"#{a}-{i}\n"+chr(39)*3+"\n = \n" + chr(39)*3

for a,b in sorted(d.items(), key=lambda x:x[0]):
	with open( f"{a}.py", "a+") as f:
		for i in range(1,b+1):
			f.write(ter(a, i)+'\n')