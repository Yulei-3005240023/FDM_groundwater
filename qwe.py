import numpy as np
import matplotlib.pyplot as plt
import hashlib
import psutil
import os

'''
a = [1, 2, 3, 4, 5]
b = [7, 7, 7, 7, 7]
c = [8, 8, 8, 8, 8]
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot()
ax.plot(a, b, linewidth=1, antialiased=True, color='red', label='jipa')
ax.plot(a, c, linewidth=1, antialiased=True, color='green')
ax.set_ylim(label='jiba')
plt.show()'''
cpu_c = psutil.cpu_count(logical=False)
path_lib = os.getcwd()
print(path_lib)
m = hashlib.md5()
m.update("1680791500.888276|1005203115".encode()) #3efa1db938bf7207cb15c68c035a929a
#resultBytes = m.digest()
resultHex = m.hexdigest()
print(resultHex)
