# PrintX

This library offers implementations of multiple modeling strategies to make a 3D model to be better printed. While the 3D printing technology has become increasingly popula3 in recent years, it still suffers from some critical limitations: small build volume, expensive printing material and long printing time. This library currently provides approaches for hollowing the 3D model into a shell, partitioning an object based on surface segmentation, etc. It is still under construction and more features will be added including strength analysis and optimization of 3D models, interlocking structure designs and 3D packing. 

The library was only tested on Windows 7 & 8. However, it should also be possible to use it on Linux and Max OS. It's easy to use since all external dependencies are included. Demos are provided to show how to use the library.

### Author
Miaojun Yao (yao.210@osu.edu)

### Version
1.0.0

### License
MIT

### Changes
* Implemented level-set-based partitioning of 3D models, as well as shell generation
* Added Demos: manual segment of a surface mesh, shell construction, partitioning

### Features
* Manual segmentation of surface mesh
* Shell generation
* Partitioning of 3D models/shells
* To be added...

### References
* Osher, S., & Fedkiw, R. (2006). Level set methods and dynamic implicit surfaces (Vol. 153). Springer Science & Business Media.
* Losasso, F., Shinar, T., Selle, A., & Fedkiw, R. (2006, July). Multiple interacting liquids. In ACM Transactions on Graphics (TOG) (Vol. 25, No. 3, pp. 812-819). ACM.
* Yao, M., Chen, Z., Luo, L., Wang, R., & Wang, H. (2015). Level-set-based partitioning and packing optimization of a printable model. ACM Transactions on Graphics (TOG), 34(6), 214.
