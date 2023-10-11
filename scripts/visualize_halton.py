import numpy as np
import matplotlib.pyplot as plt

new_u = []
new_v = []
for i in range(1, 200):
    x = 0.0
    y = 0.0
    fx = 0.5
    fy = 1 / 3
    ix = i
    while ix > 0:
        x += fx * (ix % 2)
        fx /= 2
        ix /= 2
    print(x)
    iy = i
    while iy > 0:
        y += fy * (iy % 3)
        fy /= 3
        iy /= 3
    print(y)
    r = np.sqrt(x) * 0.5
    x = r * np.sin(y * np.pi * 2)
    y = r * np.cos(y * np.pi * 2)
    new_u.append(0.5 + x * 0.020 + y * 0.025)
    new_v.append(0.5 + x * (-0.015) + y * 0.021)

print(new_u)
print(new_v)
plt.scatter(new_u, new_v, s=0.1)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.show()
