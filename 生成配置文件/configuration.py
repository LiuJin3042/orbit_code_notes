"""
修改平衡文件eqs.f
"""
# numeric: 平衡的类型,numeric = 1, 数值平衡
# numeric = 0, 解析平衡, 需要提供设备参数, 包括rmaj和a
numeric = 0
# 大径, 单位: cm
rmaj = 165
# 小径, 单位: cm
a = 80

# 波纹,krip
# krip = 0-无波纹; krip = 1-TFTR
# krip = 2-Tore Supra; krip = 3-ITER
# krip = 4-NSTX; krip = 5-Ignitor
# krip = 6-EAST
krip = 6

# 安全因子, q = q0 + qr2*r^2 + qr3*r^3, r是归一化的r
# q0: r = 0时的q值
# qed: r = 1时的q值
# qrx: r = rx时的q值为qrx. 需要设置qrx.
# 通过以上三组值解出qr2, qr3.
q0 = 1.3
qed = 4.9
rx = 0.5
qrx = 2

"""
修改扰动文件perturb.f
"""
# modes: 模数, 如果modes = 2, 可以理解为多个模叠加
# 多个模要将值写在数组里
modes = 1
mmod = [2]
nmod = [1]
amp = [5e-4]
omegv = []