from vakyume_2025_222 import export_unfinished as eu
from make_vakyume import Solver
o = eu()
S = Solver()

class Oowoo:
    def __init__(s):s.despapaye = []
    def next_end(se,x):
        '''proxima clase'''
        i = 0
        for i,l in enumerate(se.despapaye[x+1:]):
            if 'class' in l:
                return i+x + 2
            if 'import' in l: #yeah this is due to metaprogramming, pushing the script onto the end of the full-class dump
                return i+x + 2 
        return len(se.despapaye)-1

    def fugly(se,v,v_eqn):
        print(v,v_eqn)
        with open("vakyume_2025_222.py") as f:
            se.despapaye = list(f.readlines())
            aaa = f'def {v_eqn.replace("-","_")}__{v}'
            es_clase = 0
            for il in enumerate(se.despapaye):
                i,l =il
                if 'class' in l:
                    es_clase = i
                if aaa in l:
                    print(l)
                    print(es_clase,se.next_end(es_clase))
                    break

for a in o.items():
    ooo = Oowoo()
    if type(a[1])!=bool:
        # print(a[0]+'\n\t', '\n\t'.join(map(str,a[1].items())))
        print(a[0])
        for aa in sorted(a[1].items()):
            ooo.fugly(aa[0], aa[1].split(' ')[0])