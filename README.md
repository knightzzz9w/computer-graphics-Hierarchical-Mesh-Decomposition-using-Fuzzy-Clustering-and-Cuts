# computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts

C++ implements of 

[Hierarchical-Mesh-Decomposition]: https://dl.acm.org/doi/10.1145/882262.882369



## Dependence

CMake

Eigen

gtest



## Run

```bash
mkdir build
```

```bash
cd build
```

```bash
cmake ..
```

```bash
make
```

```bash
cd ..
```

```bash
./segmentNet
```



## Parameters



fuzzy-threshold: the difference between the normalized top 2 probability ; if not exceeds, then the vertice belongs to the fuzzy part



stop_dis_ratio: the max distance between the representatives divide the max distance of the origin mesh; if exceeds, the segmentation stops



stop_ang_ratio: the difference between the max dihedral angle and the min dihedral divide the max dihedral angle. if not exceeds, the segmentation stops





## Examples





- stop_dis_ratio = 0.03 ; stop_ang_ratio = 0.05

![2024-01-07 05-25-03 的屏幕截图](./examples/2024-01-07 05-25-03 的屏幕截图.png)

![2024-01-07 05-25-15 的屏幕截图](./examples/2024-01-07 05-25-15 的屏幕截图.png)

![2024-01-07 05-26-57 的屏幕截图](./examples/2024-01-07 05-26-57 的屏幕截图.png)

![2024-01-07 05-27-03 的屏幕截图](./examples/2024-01-07 05-27-03 的屏幕截图.png)

![2024-01-07 05-28-26 的屏幕截图](./examples/2024-01-07 05-28-26 的屏幕截图.png)

![2024-01-07 05-28-43 的屏幕截图](./examples/2024-01-07 05-28-43 的屏幕截图.png)



- stop_dis_ratio = 0.2 ; stop_ang_ratio = 0.3



![2024-01-07 05-25-58 的屏幕截图](./examples/2024-01-07 05-25-58 的屏幕截图.png)

![2024-01-07 05-26-04 的屏幕截图](./examples/2024-01-07 05-26-04 的屏幕截图.png)

![2024-01-07 05-27-14 的屏幕截图](./examples/2024-01-07 05-27-14 的屏幕截图.png)

![2024-01-07 05-27-49 的屏幕截图](./examples/2024-01-07 05-27-49 的屏幕截图.png)

![2024-01-07 05-29-02 的屏幕截图](./examples/2024-01-07 05-29-02 的屏幕截图.png)

![2024-01-07 05-29-18 的屏幕截图](./examples/2024-01-07 05-29-18 的屏幕截图.png)





- stop_dis_ratio = 0.4 ; stop_ang_ratio = 0.5

  

![2024-01-07 05-26-30 的屏幕截图](./examples/2024-01-07 05-26-30 的屏幕截图.png)

![2024-01-07 05-26-37 的屏幕截图](./examples/2024-01-07 05-26-37 的屏幕截图.png)

![2024-01-07 05-28-04 的屏幕截图](./examples/2024-01-07 05-28-04 的屏幕截图.png)

![2024-01-07 05-28-13 的屏幕截图](./examples/2024-01-07 05-28-13 的屏幕截图.png)

![2024-01-07 05-29-31 的屏幕截图](./examples/2024-01-07 05-29-31 的屏幕截图.png)

![2024-01-07 05-29-37 的屏幕截图](./examples/2024-01-07 05-29-37 的屏幕截图.png)
