# 减小get_n_compute_maxes_p1.f90内存使用

## 内存使用改动
- `dist_tab (n_centers_integrals, n_max_batch_size)`
  - 去掉
- `dir_tab (3,n_centers_integrals, n_max_batch_size)`
  - 去掉
- `dist_tab_sq (n_centers_integrals, n_max_batch_size)`
  - 减小维度变为 `dist_tab_sq(n_centers_integrals)`
  - 改为allocate分配内存  

## 精度对比  

||gidx - gidx'|gidx - lidx|gidx - lidx'|
|:----|:----:|:----:|:----:| 
|H2|0|1.33E-14|1.07E-14|
|H2O| 0|2.345E-12|1.927E-12|
|RBDmini| 在跑|在跑|在跑|
