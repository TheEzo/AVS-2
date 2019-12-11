/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  Tomas Willaschek <xwilla00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    10/12/2019
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.

    // 1. Compute total number of cubes in the grid.
    // todo delit na 8 podle mGridSize/2
    // todo spocitat povrch krychle > mIsoLevel+sqrt(3)/2*mGridSize
    size_t totalCubesCount = mGridSize*mGridSize*mGridSize;

    unsigned totalTriangles = count(field, mGridSize); //todo offset?

//    // 2. Loop over each coordinate in the 3D grid.
//    for(size_t i = 0; i < totalCubesCount; ++i)
//    {
//        // 3. Compute 3D position in the grid.
//        Vec3_t<float> cubeOffset( i % mGridSize,
//                                  (i / mGridSize) % mGridSize,
//                                  i / (mGridSize*mGridSize));
//
//        // 4. Evaluate "Marching Cube" at given position in the grid and
//        //    store the number of triangles generated.
//        totalTriangles += buildCube(cubeOffset, field);
//    }

    // 5. Return total number of triangles generated.
    return totalTriangles;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    for(unsigned i = 0; i < count; ++i)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        value = std::min(value, distanceSquared);
    }
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{

}

unsigned TreeMeshBuilder::count(const ParametricScalarField &field, size_t side_len) {

    return 0;
}