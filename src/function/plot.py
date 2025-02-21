import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 原始数据
Depth = np.array([5.6, 8.1, 11.3, 15.5, 20.6, 26, 31.2, 34, 37.6, 40])
amoA = np.array([1.08e7, 1.41e6, 2.51e6, 2.93e6, 2.52e6, 7.14e5, 3.63e6, 1.01e7, 4.39e6, 5.81e5])
hzo = np.array([2.70e6, 4.54e6, 4.43e6, 5.86e6, 6.51e6, 7.87e6, 4.62e6, 9.27e4, 4.45e5, 3.69e5])



# 定义指数衰减模型
def exponential_decay(z, A, B, C):
    return A * np.exp(-B * z) + C

# 拟合参数
params, _ = curve_fit(exponential_decay, Depth, amoA, p0=[1.08e7, 0.15, 5.81e5])
A, B, C = params

# 拟合结果
print(f"Fitted parameters: A = {A:.2e}, B = {B:.2f}, C = {C:.2e}")

# 生成拟合曲线
z_fit = np.linspace(5, 40, 100)
amoA_fit = exponential_decay(z_fit, A, B, C)

# 绘图
plt.figure(figsize=(10, 6))
plt.scatter(Depth, amoA, color='blue', label='Original Data')
plt.plot(z_fit, amoA_fit, color='red', label=f'Fitted Curve: $1.08e7 * e^(-0.15z) + 5.81e5$')
plt.xlabel('Depth (m)')
plt.ylabel('amoA Gene Abundance')
plt.title('amoA Gene Abundance vs Depth (Exponential Fit)')
plt.legend()
plt.grid(True)
plt.show()