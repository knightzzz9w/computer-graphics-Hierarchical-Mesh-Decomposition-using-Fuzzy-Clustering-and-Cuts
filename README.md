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

![2024-01-07 05-28-43 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/bd73939e-da13-4345-8bfe-ba478c9479a6)
![2024-01-07 05-28-26 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/27ce6284-6a14-46e8-9f20-27cb9ada406c)
![2024-01-07 05-27-03 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/6c0e12e9-85e0-4f02-a28c-c38a3f668665)
![2024-01-07 05-26-57 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/356686e4-58ae-4c2e-bce5-e933a5f17773)
![2024-01-07 05-25-15 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/da7667cc-9c81-4b58-914d-066350c0d5e1)
![2024-01-07 05-25-03 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/9b515e4c-889e-4850-b99c-70ce9c192cdb)



- stop_dis_ratio = 0.2 ; stop_ang_ratio = 0.3

![2024-01-07 05-29-18 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/91586ed6-b3cc-4283-b63c-7096abde2c95)
![2024-01-07 05-29-02 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/766546a7-9d16-4340-ae0b-a88d82db499e)
![2024-01-07 05-27-49 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/106b1da0-77e7-4cc1-bddd-2bb7a828e87a)
![2024-01-07 05-27-14 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/4f75319f-feb3-4603-b950-fe74e92f6934)
![2024-01-07 05-26-04 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/a4852266-ae36-4c2a-9134-92794a799ff4)
![2024-01-07 05-25-58 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/3a200663-961b-4333-adf9-ad3c863694c3)





- stop_dis_ratio = 0.4 ; stop_ang_ratio = 0.5

  
![2024-01-07 05-29-37 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/e619ad40-71a5-4a09-adf8-ce4c7f557ec3)
![2024-01-07 05-29-31 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/7f5899bb-7050-4ba9-b097-935eca8a5299)
![2024-01-07 05-28-13 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/60c4c8fe-632a-43ee-b050-c1ff1c1a4cc6)
![2024-01-07 05-28-04 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/d8e84942-439f-4132-8735-f93525cea90e)
![2024-01-07 05-26-37 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/8e94772e-1936-4837-b018-1c2e36014eb8)
![2024-01-07 05-26-30 的屏幕截图](https://github.com/knightzzz9w/computer-graphics-Hierarchical-Mesh-Decomposition-using-Fuzzy-Clustering-and-Cuts/assets/100104594/05be1e1f-e32c-4816-99f0-c7c3f43941b7)

