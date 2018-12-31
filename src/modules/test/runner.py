from subprocess import call

with open('results.csv', 'w') as file:
	file.write("N,R,S,avg_time,std_dev_time")

possible_N = [3, 6, 12, 24, 48, 96, 192]

for N in possible_N:
	for r in range(3, N+1):
		for s in range(3, N+1):
			print("Running ", N, r, s)
			call(['./tests', '{N}'.format(N=N), '{R}'.format(R=r), '{S}'.format(S=s)])
