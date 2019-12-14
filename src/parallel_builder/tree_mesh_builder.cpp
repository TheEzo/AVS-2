#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
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
    #pragma omp critical
    mTriangles.push_back(triangle);
}

unsigned TreeMeshBuilder::count(const ParametricScalarField &field, size_t side_len, const Vec3_t<float> start_pos) {
    if(side_len == 1){
        return buildCube(start_pos, field);
    }
    else {
        size_t side = side_len / 2; // a-1; a
        CubeCornerVerts_t cubeCorners;
        Vec3_t<float> to_count(side, side, side);
        transformCubeVertices(to_count, sc_vertexNormPos, cubeCorners);
        if (!continueCount(cubeCorners[0], field, side))
            return 0;

        unsigned t1 = 0;
        unsigned t2 = 0;
        unsigned t3 = 0;
        unsigned t4 = 0;
        unsigned t5 = 0;
        unsigned t6 = 0;
        unsigned t7 = 0;
        unsigned t8 = 0;
        #pragma omp parallel
        {
            #pragma omp master
            {
                #pragma omp task //private(t1)
                {
                    Vec3_t<float> start(start_pos.x, start_pos.y, start_pos.z);
                    t1 = count(field, side, start);
                }
                #pragma omp task //private(t2)
                {
                    Vec3_t<float> start(start_pos.x + side, start_pos.y, start_pos.z);
                    t2 = count(field, side, start);
                }
                #pragma omp task //private(t3)
                {
                    Vec3_t<float> start(start_pos.x, start_pos.y + side, start_pos.z);
                    t3 = count(field, side, start);
                }
                #pragma omp task //private(t4)
                {
                    Vec3_t<float> start(start_pos.x, start_pos.y, start_pos.z + side);
                    t4 = count(field, side, start);
                }
                #pragma omp task //private(t5)
                {
                    Vec3_t<float> start(start_pos.x + side, start_pos.y + side, start_pos.z);
                    t5 = count(field, side, start);
                }
                #pragma omp task //private(t6)
                {
                    Vec3_t<float> start(start_pos.x, start_pos.y + side, start_pos.z + side);
                    t6 = count(field, side, start);
                }
                #pragma omp task //private(t7)
                {
                    Vec3_t<float> start(start_pos.x + side, start_pos.y, start_pos.z + side);
                    t7 = count(field, side, start);
                }
                #pragma omp task //private(t8)
                {
                    Vec3_t<float> start(start_pos.x + side, start_pos.y + side, start_pos.z + side);
                    t8 = count(field, side, start);
                }
            }
        }
        #pragma omp taskwait
        return t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8;
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

#pragma clang diagnostic pop