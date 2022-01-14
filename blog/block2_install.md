# block2安装
## 准备工作
下载block2
```
git clone https://github.com/block-hczhai/block2-preview
```
安装pybind11
```
pip bybind11
```
加载CMake 3.7
```
module add cmake/3.7.0
```
加载gcc 8.3
```
module add gcc/8.3.0-wzm
```
## 开始安装
可以在 <block2路径>/setup.py 的cmake_args变量中配置cmake参数，具体参数可以在README.md中查看。我增加了开启mpi的参数。

执行命令编译并创建whl包
```
python3 setup.py bdist_wheel
```
创建的whl包默认在 <block2路径>/dist 文件夹下，使用pip安装whl包
```
pip install block2-0.1.10-cp37-cp37m-linux_x86_64.whl
```
安装过程中将自动检查依赖包并下载。
