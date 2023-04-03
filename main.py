import numpy as np
from scipy.stats import bernoulli
from scipy.stats import entropy
import matplotlib.pyplot as plt
import random
import collections


def getEntropy(arr):
    amounts = collections.Counter([value for value in arr])
    distributions = [x / sum(amounts.values()) for x in amounts.values()]
    return entropy(distributions, base=2)


# Bernoulli
def bernoulli_generate(p, n):
    rand_arr = [random.random() for _ in range(n)]
    return [1 if rand < p else 0 for rand in rand_arr]


# 12 % 3 = 0 - Binomial
def binomial_generate(p, n):
    result = []
    for i in range(n):
        x = 0
        for j in range(100):
            if random.random() < p:
                x += 1
        result.append(x)
    return result


N = 23
p = 0.5 + (N - 12) / 1000
start = 10
end = 1000
step = 10
selectionLen = 1000

# default entropy
entropy_bern = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
entropy_binom = 1/2 * np.log2(2*np.pi*np.e*100*p*(1-p))
print(entropy_bern)
print(entropy_binom)

dispersion_BRV_bern = []
dispersion_BRV_binom = []
dispersion_built_in_bern = []
dispersion_built_in_binom = []

for _len in range(start, end, step):
    deviation_BRV_bern = []
    deviation_BRV_binom = []
    deviation_built_in_bern = []
    deviation_built_in_binom = []

    for n in range(selectionLen):
        BRV_bern = bernoulli_generate(p, _len)
        BRV_binom = binomial_generate(p, _len)
        built_in_bern = bernoulli.rvs(p, size=_len)
        built_in_binom = np.random.binomial(_len, p, _len)

        deviation_BRV_bern.append(entropy_bern - getEntropy(BRV_bern))
        deviation_BRV_binom.append(entropy_binom - getEntropy(BRV_binom))
        deviation_built_in_bern.append(entropy_bern - getEntropy(built_in_bern))
        deviation_built_in_binom.append(entropy_binom - getEntropy(built_in_binom))

    dispersion_BRV_bern.append(np.var(deviation_BRV_bern))
    dispersion_BRV_binom.append(np.var(deviation_BRV_binom))
    dispersion_built_in_bern.append(np.var(deviation_built_in_bern))
    dispersion_built_in_binom.append(np.var(deviation_built_in_binom))

plt.subplot(2, 2, 1)
plt.plot(range(start, end, step), dispersion_BRV_bern)
plt.title('Bernoulli BRV distribution')
plt.show()

plt.subplot(2, 2, 2)
plt.plot(range(start, end, step), dispersion_built_in_bern)
plt.title('Bernoulli built-in distribution')
plt.show()

plt.subplot(2, 2, 3)
plt.plot(range(start, end, step), dispersion_BRV_binom)
plt.title('Binomial BRV distribution')
plt.show()

plt.subplot(2, 2, 4)
plt.plot(range(start, end, step), dispersion_built_in_binom)
plt.title('Binomial built-in distribution')
plt.show()