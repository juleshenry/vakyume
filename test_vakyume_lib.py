import vakyume_lib
if __name__=='__main__':
	for o in dir(vakyume_lib):
		stro = str(o)
		if stro[0].isalpha():
			print('sublib',stro)
			sublib = getattr(vakyume_lib,stro)
			for a in dir(sublib):
				astro = str(a)
				if stro[0].isalpha():
					x = vars(vakyume_lib)[stro]
					print(getattr(x,astro))

					# astro = str(a)
					# sublibmethod = getattr(o, astro)
					# print(sublibmethod)
		break
