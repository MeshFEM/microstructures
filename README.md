# Microstructures

This repository contains a collection of routines for the optimization and analysis of homogenized material properties of graph-based microstructures. Microstructures are described by a graph, with parameters such as vertex positions, vertex radii, and blending factor. They can be inflated into a periodic triangle mesh (2D) or tetrahedral mesh (3D).

### Coverage Results in 2D

See this [coverage plot](http://julianpanetta.com/2d_isosurface_hull_autocover_results/flipper.html).

### Dependencies

The following dependencies are required, and downloaded automatically by CMake (expect boost):
- [MeshFEM](https://github.com/MeshFEM/MeshFEM)
- [TBB](https://github.com/wjakob/tbb)
- [CLI11](https://github.com/wjakob/tbb)
- [Nanoflann](https://github.com/jlblancoc/nanoflann)
- [Quickhull](https://github.com/akuukka/quickhull.git)
- [Catch2](https://github.com/catchorg/Catch2)
- [Boost](https://www.boost.org/) (installed separately)

The following dependencies are optional:
- [CGAL](https://github.com/CGAL/cgal)
- [Ceres](https://github.com/ceres-solver/ceres-solver)
- [Dlib](https://github.com/davisking/dlib)
- [Knitro](https://www.artelys.com/solvers/knitro/)
- [Libigl](https://github.com/libigl/libigl.git)
- [Cotire](https://github.com/sakra/cotire.git)
- [Sanitizers](https://github.com/arsenm/sanitizers-cmake.git)

### Compilation

```
mkdir build && cd build
cmake ..
make -j 8
```

### Content

#### Libraries

This repository is organized into separate libraries, located under `src/lib`:
- **isosurface_inflator**: Generates a triangle/tet mesh from graph labeled with vertex thicknesses.
- **inflators**: Wrappers around various inflators (including the isosurface_inflator) providing a unified API.
- **optimizers**: Wrappers around various numerical optimization packages for minimizing a constrained non-linear problem.
- **pattern_optimization**: Optimizes the parameters of a particular microstructure to achieve a target material property.

#### Binaries

The functionalities are fragmented into a various executables that can be used to enumerate different graph topology, inflate a given microstructure, optimize parameters for a target elastic tensor or minimize worst-case stress.
Please refer to individual binary usages for additional details.

### Usage

Documentation on how to run worst-case stress optimization is available here:
- [Worst-Case Stress Optimization](docs/usage.md)
- [Common Optimization Problems](docs/notes.md)
- [FAQ](docs/faq.md)

### Citation

The code in this repository was used in the following papers:

<details><summary>Panetta, J., Zhou, Q., Malomo, L., Pietroni, N., Cignoni, P. and Zorin, D., 2015. Elastic textures for additive fabrication. ACM Transactions on Graphics (TOG), 34(4), pp.1-12.</summary>

```bibtex
@article{Panetta2015,
    author = {Panetta, Julian and Zhou, Qingnan and Malomo, Luigi and Pietroni, Nico and Cignoni, Paolo and Zorin, Denis},
    title = {Elastic Textures for Additive Fabrication},
    journal = {ACM Trans. Graph.},
    issue_date = {August 2015},
    volume = {34},
    number = {4},
    month = jul,
    year = {2015},
    issn = {0730-0301},
    pages = {135:1--135:12},
    articleno = {135},
    numpages = {12},
    disabled_url = {http://doi_disabled_disabled.acm.org/10.1145/2766937},
    doi_disabled_disabled = {10.1145/2766937},
    acmid = {2766937},
    publisher = {ACM},
    address = {New York, NY, USA},
    keywords = {additive fabrication, deformable objects, goal-based material design, homogenization, microstructures, shape optimization},
}
```
</details>

<details><summary>Panetta, J., Rahimian, A. and Zorin, D., 2017. Worst-case stress relief for microstructures. ACM Transactions on Graphics (TOG), 36(4), pp.1-16.</summary>

```bibtex
@article{Panetta2015,
    author = {Panetta, Julian and Zhou, Qingnan and Malomo, Luigi and Pietroni, Nico and Cignoni, Paolo and Zorin, Denis},
    title = {Elastic Textures for Additive Fabrication},
    journal = {ACM Trans. Graph.},
    issue_date = {August 2015},
    volume = {34},
    number = {4},
    month = jul,
    year = {2015},
    issn = {0730-0301},
    pages = {135:1--135:12},
    articleno = {135},
    numpages = {12},
    disabled_url = {http://doi_disabled_disabled.acm.org/10.1145/2766937},
    doi_disabled_disabled = {10.1145/2766937},
    acmid = {2766937},
    publisher = {ACM},
    address = {New York, NY, USA},
    keywords = {additive fabrication, deformable objects, goal-based material design, homogenization, microstructures, shape optimization},
}
```
</details>

<details><summary>Tozoni, D.C., Dumas, J., Jiang, Z., Panetta, J., Panozzo, D. and Zorin, D., 2020. A low-parametric rhombic microstructure family for irregular lattices. ACM Transactions on Graphics (TOG), 39(4), pp.101-1.</summary>

```bibtex
@article{Tozoni2020,
  doi = {10.1145/3386569.3392451},
  url = {https://doi.org/10.1145/3386569.3392451},
  year = {2020},
  month = jul,
  publisher = {Association for Computing Machinery ({ACM})},
  volume = {39},
  number = {4},
  author = {Davi Colli Tozoni and J{\'{e}}r{\'{e}}mie Dumas and Zhongshi Jiang and Julian Panetta and Daniele Panozzo and Denis Zorin},
  title = {A low-parametric rhombic microstructure family for irregular lattices},
  journal = {{ACM} Transactions on Graphics}
}
```
</details>


### Licensing

This project is licensed under the MIT License. See [LICENSE](LICENSE) for more information.

**Exception**: Some files in this repository depend on CGAL, which is licensed under GPL. Those parts are disabled by default, but can be enabled with the preprocessor macro `MICRO_WITH_CGAL`, by turning on the corresponding CMake option. The parts of this code that depend on CGAL are licensed under GPL. See [LICENSE.GPL](LICENSE.GPL) for more information.

### Acknowledgment

The main author of this library is Julian Panetta, with contributions from Jérémie Dumas, Davi Colli Tozoni, Luigi Malomo, Qingnan (James) Zhou.
