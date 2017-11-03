def prettyDir(obj):
  return {"members":",".join(getattr(obj,'__members__',[])),
      "methods":"(),".join(getattr(obj,'__methods__',[]))+"()"}

members=["vertices","edges","faces","cells","nodes","elements"]
d={}
for mb in members:
  objArray=getattr(part1,mb)
  d[mb]=prettyDir(objArray)
  obj=objArray[0]
  d[mb[:-1]]=prettyDir(obj)
for k,v in d.items():
  print k,":",v