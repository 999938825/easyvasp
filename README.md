# easyvasp
高通量建立vasp模型和批量分析模型性质
(High throughput VASP model and batch analysis model properties)
这个项目重新封装了pymatgen的部分功能,以便非专业的研究者快速使用高通量建模
(This project repackaged some of Pymatgen's capabilities to allow non-specialist researchers to quickly use high-throughput modeling)

快速开始

首先安装c++14 否则pymatgen安装容易报错
安装Anaconda3 作为python环境
安装pymatgen及其依赖包，建议使用国内源
pip install pymatgen==2018.8.13 -i https://pypi.mirrors.ustc.edu.cn/simple/
安装scikit-learn，pandas 方法同上
安装easyvasp
1.  下载赝势文件夹，并复制到 easyvasp\easyvasp\Pseudopotentials 文件夹，赝势文件夹的内部结构类似于pbe\Co 这样
1.	打开cmd
2.	cd /d “easyvasp安装目录”
3.	Python setup.py install 安装easyvasp
4.	任意安装一个python IDE 现在你就可以编写程序来建模了！（推荐pycharmm）
教程请参考easyvasp教程.docx
