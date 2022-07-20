ITKFpfh
=================================

[<img src="https://github.com/PranjalSahu/Fpfh.git/actions/workflows/build-test-package.yml/badge.svg">](https://github.com/PranjalSahu/Fpfh/actions)
[<img src="https://img.shields.io/pypi/v/itk-fpfh.svg">](https://pypi.python.org/pypi/itk-fpfh)


Overview
--------

Module to calculate FPFH feature for a pointset.
Sample Usage is shown below:

```
# normal_pointset is ITK Pointset which contains normal vector for each point
# pointset is ITK Pointset which contains the input points for which feature needs to be calculated

# normal_np is numpy array of shape [Nx3]
# fpfh_feature is numpy array of shape [Nx33]

normal_pointset.SetPoints(itk.vector_container_from_array(normal_np.flatten()))
fpfh = itk.Fpfh.PointFeature.MF3MF3.New()
fpfh.ComputeFPFHFeature(pointset, normal_pointset, 25, 100)
fpfh_feature = fpfh.GetFpfhFeature()
fpfh_feature = itk.array_from_vector_container(fpfh_feature)
```

One can obtain the normals using the following code:
```
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
