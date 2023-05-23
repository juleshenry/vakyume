import vakyume_lib
if __name__=='__main__':
	for o in dir(vakyume_lib):
		stro = str(o)
		if stro[0].isalpha():
			sublib = getattr(vakyume_lib,stro)
			for a in dir(sublib):
				astro = str(a)
				sublibmethod = getattr(str(astro))
		
