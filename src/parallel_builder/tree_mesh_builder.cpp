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
    Vec3_t<float> start_pos(0, 0, 0);
    unsigned totalTriangles = count(field, mGridSize, start_pos);

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
        mTriangles.push_back(triangle);
}

unsigned TreeMeshBuilder::count(const ParametricScalarField &field, size_t side_len, const Vec3_t<float> start_pos) {
    if(side_len == 1){
        return buildCube(start_pos, field);
    }
    else{
        size_t side = side_len/2; // a-1; a
        CubeCornerVerts_t cubeCorners;
        Vec3_t<float> to_count(side, side, side);
        transformCubeVertices(to_count, sc_vertexNormPos, cubeCorners);
        if(!continueCount(cubeCorners[0], field, side))
            return 0;

        unsigned total_triangles = 0;
        {
            Vec3_t<float> start(start_pos.x, start_pos.y, start_pos.z);
            total_triangles += count(field, side, start);
        }

        {
            Vec3_t<float> start(start_pos.x + side, start_pos.y, start_pos.z);
            total_triangles += count(field, side, start);
        }

        {
            Vec3_t<float> start(start_pos.x, start_pos.y + side, start_pos.z);
            total_triangles += count(field, side, start);
        }

        {
            Vec3_t<float> start(start_pos.x, start_pos.y, start_pos.z + side);
            total_triangles += count(field, side, start);
        }

        {
            Vec3_t<float> start(start_pos.x + side, start_pos.y + side, start_pos.z);
            total_triangles += count(field, side, start);
        }

        {
            Vec3_t<float> start(start_pos.x, start_pos.y + side, start_pos.z + side);
            total_triangles += count(field, side, start);
        }

        {
            Vec3_t<float> start(start_pos.x + side, start_pos.y, start_pos.z + side);
            total_triangles += count(field, side, start);
        }

        {
            Vec3_t<float> start(start_pos.x + side, start_pos.y + side, start_pos.z + side);
            total_triangles += count(field, side, start);
        }
        return total_triangles;
    }

}

bool TreeMeshBuilder::continueCount(const Vec3_t<float> pos, const ParametricScalarField &field, size_t side){
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
    return sqrt(value) <= mIsoLevel+sqrt(3)/2*side;
}
