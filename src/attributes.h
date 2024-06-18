#pragma once

#include <Eigen/Core>

class VertexAttributes
{
public:
    VertexAttributes(double x = 0, double y = 0, double z = 0, double w = 1)
    {
        position << x, y, z, w;
        xyz << x, y, z;
    }

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes &a,
        const VertexAttributes &b,
        const VertexAttributes &c,
        const double alpha,
        const double beta,
        const double gamma)
    {
        VertexAttributes r;
        r.position = alpha * a.position + beta * b.position + gamma * c.position;
        r.normal = alpha * a.normal + beta * b.normal + gamma * c.normal;
        r.color = alpha * a.color + beta * b.color + gamma * c.color;
        r.xyz = alpha * a.xyz + beta * b.xyz + gamma * c.xyz;
        r.triangles_touching = 0;
        return r;
    }

    Eigen::Vector4d position;
    Eigen::Vector4d color;
    Eigen::Vector3d normal;
    Eigen::Vector3d xyz;
    int triangles_touching;
};

class FragmentAttributes
{
public:
    FragmentAttributes(double r = 0, double g = 0, double b = 0, double a = 1)
    {
        color << r, g, b, a;
    }

    Eigen::Vector4d color;
    Eigen::Vector4d position;

};

class FrameBufferAttributes
{
public:
    FrameBufferAttributes(double r = 0, double g = 0, double b = 0, double a = 1)
    {
        color << r, g, b, a;
        depth = -1;
    }

    Eigen::Matrix<double, 4, 1> color;
    double depth;

};

class UniformAttributes
{
public:
    Eigen::Matrix4d view;
    Eigen::Vector4d color;
    Eigen::Matrix4d projective;
    Eigen::Matrix4d cam;
    Eigen::Matrix4d orthographic;
    Eigen::Matrix4d perspective;
    
};