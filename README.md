# 3DKMI - Three-Dimensional Krawtchouk Moment Invariants

This is the **MATLAB** package for extracting **3D Krawtchouk moment invariants**.

It can extract both local and global features from 3D images.

The features derived from the 3D images remain stable under translation, scaling, and rotation.

The program’s dependencies include the `Symbolic Math Toolbox` and the `Image Processing Toolbox`.

## 🧪 Example Usage
```
dim = 128;
path = 'please input the path that contains 2D image slices or a single 3D image file';
arr = readImage(path, dim);
invariants = get3DKMI(arr, dim);
```
Two 3D images from McGill's 3D shape benchmark website and one example data of *Menardella exilis* in the `example_data` folder are provided for code execution.

A quick-start guide and more examples are available at `https://krawtchouk.github.io/`.

## 📝 Citation
If you use this package in your research, please cite it.
```
@article{lin20243dkmi,
  title={3DKMI: A MATLAB package to generate shape signatures from Krawtchouk moments and an application to species delimitation in planktonic foraminifera},
  author={Lin, Huahua and Zhang, Wenshu and Mulqueeney, James M and Brombacher, Anieke and Searle-Barnes, Alex and Nixon, Mark and Cai, Xiaohao and Ezard, Thomas HG},
  journal={Methods in Ecology and Evolution},
  volume={15},
  number={11},
  pages={1940--1948},
  year={2024},
  publisher={Wiley Online Library}
}
```
