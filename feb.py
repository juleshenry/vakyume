from vakyume_2025 import export_unfinished as eu
o = eu()
for a in o.items():
    if type(a[1])!=bool:print(a[0]+'\n\t', '\n\t'.join(map(str,a[1].items())))