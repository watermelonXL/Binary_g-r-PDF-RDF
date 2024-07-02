# Binary_g-r-PDF-RDF
# 这是一个关于计算2元体系RDF(radial pair distribution function g(r))或者PDF的代码。
## 文件介绍 

**PDF_PPDF.INI**      这是一个配置文件，需要计算的参数都在里面。 

**PDF_PPDF.f90**      这是Fortran程序的源码。 

**PDF_PPDF.exe**      这是一个可执行程序。

# 如何计算PDF呢？这里需要几个步骤。
## 第一步：确定需要分析的文件格式，本程序只能分析如下相似的格式，
```
ITEM: TIMESTEP
10500000
ITEM: NUMBER OF ATOMS
50460
ITEM: BOX BOUNDS pp pp pp
0.0278173 83.1007
0.0278173 83.1007
0.0287765 85.9662
ITEM: ATOMS id type x y z
17735 1 1.96208 1.46122 2.76805
6746 1 2.01231 3.23008 1.0509
32073 1 3.73142 0.90597 1.26289
```
此种格式按照Lammps中dump的custom中的进行输出，详细的输出命令介绍可参照如下链接， 

Lammps的dump输出文件格式 - 我爱吃西瓜的文章 - 知乎
https://zhuanlan.zhihu.com/p/706138884 

[Lammps的dump输出文件格式 - 我爱吃西瓜的文章 - 知乎](https://zhuanlan.zhihu.com/p/706138884)

## 第二步，将需要分析的三个文件放在同一目录下。
三个文件分别是： 

**Lammps的输出文件** 

**PDF_PPDF.exe** 

**PDF_PPDF.INI** 


## 第三步，点击**PDF_PPDF.exe**，等待程序运行完毕。
