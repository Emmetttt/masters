import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, ks_2samp
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

##Create two Gaussian distributions of random numbers, differing
##only in a standard deviation of 0.2.
dist1 = np.random.normal(0, 1, 500)
dist2 = np.random.normal(0, 0.8, 2000)

print(ks_2samp(dist1, dist2)) ##Perform a KS test

dist1.sort() ##Sort the two distributions from low to high to create
dist2.sort() ##a cumulative frequency plot

##Plotting
plt.plot(dist1, fun.cumfreq(dist1), label="Distribution 1", color='blue')
plt.plot(dist2, fun.cumfreq(dist2), label="Distribution 2", color='red')
emptyNLSS = plt.hist([], range=[-0.3,0.3], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(dist1, dist2)[0], ks_2samp(dist1, dist2)[1]), color='white')
plt.xlabel('Value')
plt.ylabel('Cumulative Frequency')
plt.legend()
plt.show()

##Empirical p-value determination
bothdist = np.append(dist1, dist2) ##Parent distribution
trials = 10000
stats = []
i=0
while i < trials:
        np.random.shuffle(bothdist)
        rand1 = np.random.randint(0,2499,size=500) ##Creates 500x1 array of random numbers
        rand2 = np.random.randint(0,2499,size=2000)
        test1 = [float(bothdist[y]) for y in rand1] ##Associates random numbers with elements
        test2 = [float(bothdist[y]) for y in rand2] ##in bothdist array
        stats.append(ks_2samp(test1, test2)[0]) ##Performs KS test, stores stat
        i+=1

stats.sort() ##Sort from low to high
j = 0
statcrit = ks_2samp(dist1, dist2)[0] ##Retrieves critical p-val
while j < len(stats): ##Loop over all stats
    if stats[j] > statcrit:
        print(((len(stats)-j)/len(stats))*100, '% of randomly generated sets have stat value > ', statcrit)
        break ##Stops loop if p=pcrit
    else:
        j+=1

plt.hist(stats, bins=100)
plt.xlabel('Stat')
plt.ylabel('Frequency')
plt.show()
