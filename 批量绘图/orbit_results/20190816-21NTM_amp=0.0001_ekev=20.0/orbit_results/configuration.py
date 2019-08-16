# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 17:55:23 2019

Read the configuration file, read the file to be changed, and write the new file

@author: LJ
"""

"""
修改平衡文件eqs.f
"""

# numeric: 平衡的类型,numeric = 1, 数值平衡
# numeric = 0, 解析平衡, 需要提供设备参数, 包括rmaj和a
numeric = 0
if numeric == 0:
    # rmaj, 大径, 半径, 单位: cm
    rmaj = 185
    # a, 小径, 直径, 单位: cm
    a = 40
    # 安全因子, q = q0 + qr2*r^2 + qr3*r^3, r是归一化的r
    # r在源代码中是以rmaj归一, 但是为了便于设置, 这里还是采用以a归一
    # q0: r = 0时的q值
    # qed: r = 1时的q值
    # qrx: r = rx时的q值为qrx. 需要设置qrx.
    # 通过以上三组值解出qr2, qr3.
    q0 = 1.3
    qed = 4.9
    rx = 0.5
    qrx = 2

# 波纹种类,krip
# krip = 0-无波纹; krip = 1-TFTR
# krip = 2-Tore Supra; krip = 3-ITER
# krip = 4-NSTX; krip = 5-Ignitor
# krip = 6-EAST
krip = 0


"""
修改扰动文件perturb.f
"""
# modes: 模数, 如果modes = 2, 可以理解为多个模叠加
# 多个模要将值写在数组里
# mmod: m值
# nmod: n值
# amp: 模幅度
# omegv: 频率, 单位为千赫兹
# dele: 能量改变的步长, 只有在omegv不为0的时候才要设置
modes = 1
harm = [1,2,3]
mmod = [2,2,2]
nmod = [1,1,1]
amp = [5e-4,5e-4,5e-4]
omegv = [0,1,2]
dele = 10
alfv = [1,1,1]
        

# 扰动模式
# a1 = 1-gaussian, a1 = 2-gaussian MHD
# a1 = 3-MHD, a1 = 4-resistive
a1 = 4

""" 
修改主程序文件orbit.F
"""
# npert: 扰动模式
# npert = 1-有磁扰动, npert = 2-有势扰动
# npert = 3-两者都有, npert = 0-无扰动
npert = 1

# nplot: 运行模式
# nplot = 1-单粒子, nplot = 2-多粒子
nplot = 3

# pdist, 粒子分布模式
# pdist = 1-shelldep, pdist = 2-sampledep
# pdist = 3-poindep,  pdist = 4-poinkdep
pdist = 3

if pdist != 2:
    # nprt: 粒子数
    nprt = 50
    
    # ntor: 程序运行时间, 粒子绕环ntor周, 程序停止
    ntor = 1500
    
    # bkg: 磁场强度, 千高斯
    bkg = 18
    
    # polo, shelldep粒子起始分布磁面
    # p1, p2, poindep粒子分布的起始结束磁面
    # pchi, 粒子起始俯仰角
    polo = 0.5
    p1 = 0.01
    p2 = 0.99
    pchi = 0
    
    # zprt, 粒子带电荷数, 单位为1个单位电荷
    # prot, 粒子质量, 单位为单个质子质量
    # ekev, 粒子能量, 单位为千电子伏
    zprt = 1
    prot = 2
    ekev = 60

# comment: 对本文件的目的说明, 可以随意修改, 会被用作新生成的文件夹名称
comment = '21NTM_amp=%s'%amp[0]































