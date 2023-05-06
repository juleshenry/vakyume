import os
import re

def see_which_notes_are_valid_Python():
	for o in os.listdir(os.getcwd()+'/chapters'):
		if o[0]=='_':continue

		with open(os.getcwd()+'/chapters/'+o) as s:
			eqn_number = ""
			for l in s.readlines():
				if (x:=re.compile('\d{1,2}-\d{1,2}').findall(l)):
					eqn_number = x[0]
				if '=' in l and ':' not in l.split('=')[0]:
					parse_eqn(l.strip().split('#')[0])


def parse_eqn(l):
	try:
		eval(l)
	except SyntaxError as se:
		to_tokens = se.text.split('=')[0]
		y = re.compile('\w[A-Za-z0-9\-\_]+').findall(to_tokens)
		for o in y:
			print(o)
	except NameError as ne:
		print(ne)


if __name__=='__main__':
	see_which_notes_are_valid_Python()