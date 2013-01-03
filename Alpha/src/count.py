import fileinput
def countUses():
	instances=dict()
	for line in fileinput.input():
		lin=line.strip("\n").split(" ")
		lin.pop(len(lin)-1)
		#print(lin)
		for i in lin:
			key=int(i)
			if key in instances:
				instances[key]+=1
			else:
				instances[key]=1
	print("Number of entries = {}".format(len(instances)))
	val=instances.items()
	vallist=[]
	for i in iter(val):
		vallist.append(i)
	vallist.sort()
	for i in vallist:
		print(i)
#	print(vallist)
	
if __name__=='__main__':
	countUses();
