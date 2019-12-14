/**
 * @file    tree_mesh_builder.h
 *
 * @author  Tomas Willaschek <xwilla00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    10/12/2019
 **/

#ifndef TREE_MESH_BUILDER_H
#define TREE_MESH_BUILDER_H

#include "base_mesh_builder.h"

class TreeMeshBuilder : public BaseMeshBuilder
{
public:
    TreeMeshBuilder(unsigned gridEdgeSize);

protected:
    unsigned marchCubes(const ParametricScalarField &field);
    unsigned count(const ParametricScalarField &field, size_t side_len, const Vec3_t<float> start_pos);
    bool continueCount(const Vec3_t<float> pos, const ParametricScalarField &field, size_t side);
    float evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field);
    void emitTriangle(const Triangle_t &triangle);
    const Triangle_t *getTrianglesArray() const { return mTriangles.data(); }

    std::vector<Triangle_t> mTriangles; ///< Temporary array of triangles
};

#endif // TREE_MESH_BUILDER_H
