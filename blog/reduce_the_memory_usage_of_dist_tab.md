# 减少dist_tab,dir_tab和dist_tab_sq的内存使用

### 原始代码(c语言描述)

位置：src/DFPT_dielectric/integrate_first_order_rho_dielectric.f90和src/DFPT_dielectric/integrate_first_order_H_dielectric.f90

```c
double dist_tab_sq[n_max_batch_size][n_centers_integrals];
double dist_tab[n_max_batch_size][n_centers_integrals];
double dir_tab[n_max_batch_size][n_centers_integrals][3];

for(int i_index=0;i_index<bach_size;i_index++)
{
    //每次循环计算一行dist_tab_sq和dir_tab
    tab_atom_centered_coords_p0(dist_tab_sq[i_index],dir_tab[i_index]);
}
for(int i_index=0;i_index<batch_size;i_index++)
{
    //每次循环使用一行dist_tab_sq,dir_tab和dist_tab
    functions(dist_tab_sq[i_index],dist_tab[i_index],dir_tab[i_index]);
}
```

### 修改后的代码

可以看到这三个数组是一行一行计算得到也是一行一行使用的，每一行之间也没有联系，所以可以循环合并，计算了一行便可以使用一行，而不用把每一行都计算出来

```c
double dist_tab_sq[n_centers_integrals];
double dist_tab[n_centers_integrals];
double dir_tab[n_centers_integrals][3];

for(int i_index=0;i_index<batch_size;i_index++)
{
    tab_atom_centered_coords_p0(dist_tab_sq,dir_tab);
    functions(dist_tab_sq,dist_tab,dir_tab);
}
```

