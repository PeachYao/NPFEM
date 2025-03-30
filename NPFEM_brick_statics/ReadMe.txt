主要参考来源
1.  https://github.com/HaoguangYang/OpenSTAP/blob/master/src/stap.f90
2. https://github.com/PeachYao/STAP90


主要改编位于 brick.f90，改写了刚度矩阵S，添加了等效力向量FK
主程序中的形变U被替换为了deformed position。