ITKFpfh
=================================

[<img src="https://github.com/PranjalSahu/Fpfh/actions/workflows/build-test-package.yml/badge.svg">](https://github.com/PranjalSahu/Fpfh/actions)
[<img src="https://img.shields.io/pypi/v/itk-fpfh.svg">](https://pypi.python.org/pypi/itk-fpfh)


Overview
--------
<p align="center">
<img src="https://user-images.githubusercontent.com/1044135/180022085-0e230afa-a5df-4afe-86b5-29543a2e1823.png" width="400" height="400" />
</p>

Module to calculate FPFH feature for a pointset.
Sample Usage is shown below:

```python
# normal_pointset is ITK Pointset which contains normal vector for each point
# pointset is ITK Pointset which contains the input points for which feature needs to be calculated

# normal_np is numpy array of shape [Nx3]
# fpfh_feature is numpy array of shape [33xN]
# 25 is the radius and 100 is the maximum number of neighbors

normal_pointset.SetPoints(itk.vector_container_from_array(normal_np.flatten()))
fpfh = itk.Fpfh.PointFeature.MF3MF3.New()
fpfh.ComputeFPFHFeature(pointset, normal_pointset, 25, 100)
fpfh_feature = fpfh.GetFpfhFeature()
fpfh_feature = itk.array_from_vector_container(fpfh_feature)
fpfh_feature = np.reshape(fpfh_feature, [33, pointset.GetNumberOfPoints()])
```

One can obtain the normals using the following code:
```python
def getnormals_pca(inputPoints):
    import vtk
    from vtk.util import numpy_support
    meshPoints = numpy_to_vtk_polydata(inputPoints)
    normals = vtk.vtkPCANormalEstimation()
    normals.SetSampleSize(30)
    normals.SetFlipNormals(True)
    #normals.SetNormalOrientationToPoint()
    normals.SetNormalOrientationToGraphTraversal()
    normals.SetInputData(meshPoints)
    normals.Update()
    as_numpy = numpy_support.vtk_to_numpy(normals.GetOutput().GetPointData().GetArray(0))
    return as_numpy
```
