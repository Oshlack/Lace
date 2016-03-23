f = open('output.pslx','r')
iter = 0
for line in f:
	print(iter)
	if(iter < 4): print(line)
	iter=iter+1
