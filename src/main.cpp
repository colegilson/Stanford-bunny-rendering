// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

#include <gif.h>
#include <fstream>

#include <Eigen/Geometry>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;

//Image height
const int H = 480;

//Camera settings
const double near_plane = 1.5;       //AKA focal length
const double far_plane = near_plane * 100;
const double field_of_view = 0.7854; //45 degrees
const double aspect_ratio = 1.5;
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 3);
const Vector3d camera_gaze(0, 0, -1);
const Vector3d camera_top(0, 1, 0);

//Object
const std::string data_dir = DATA_DIR;
const std::string mesh_filename(data_dir + "bunny.off");
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)

//Material for the object
const Vector3d obj_diffuse_color(0.5, 0.5, 0.5);
const Vector3d obj_specular_color(0.2, 0.2, 0.2);
const double obj_specular_exponent = 256.0;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector3d> light_colors;
//Ambient light
const Vector3d ambient_light(0.3, 0.3, 0.3);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    if (!in.good())
    {
        std::cerr << "Invalid file " << mesh_filename << std::endl;
        exit(1);
    }
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16);
}

void build_uniform(UniformAttributes &uniform)
{
    //TODO: setup uniform

    uniform.color << 1, 0, 0, 1;

    //TODO: setup camera, compute w, u, v
    uniform.view << 1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;

    if (aspect_ratio < 1)  
        uniform.view(0,0) = aspect_ratio;
    else   
        uniform.view(1,1) = aspect_ratio;
    
    Vector3d w = -(camera_gaze.normalized()); // -g/||g||
    Vector3d u = (camera_top.cross(w).normalized()); // (t x u)/ ||t x u ||
    Vector3d v = w.cross(u); // w x u

    //TODO: compute the camera transformation
    double n = near_plane;
    double f = -far_plane;
    double t = n * tan(field_of_view / 2);
    double b = -t;
    double r = t * aspect_ratio;
    double l = -r;
    Vector3d e = camera_position;

    // the following is what is shown in the textbook, I did it this way
    // instead of what is shown in the slides as taking the inverse of 
    // [u, v, w, e]
    // [0, 0, 0, 1].inverse().eval()
    // left the uniform.cam with values of -0 which upset me :)
    // cammera matrix, Mcam
    uniform.cam <<  u.x(), u.y(), u.z(), 0,
                    v.x(), v.y(), v.z(), 0,
                    w.x(), w.y(), w.z(), 0,
                    0,     0,     0,     1;
    Matrix4d cam;
    cam <<  1, 0, 0, -e.x(),
            0, 1, 0, -e.y(),
            0, 0, 1, -e.z(),
            0, 0, 0, 1;
    uniform.cam = (uniform.cam * cam);
    //cout << uniform.cam << endl;

    // orthogonal matrix, Morth
    uniform.orthographic << 2 / (r-l), 0,         0,            -(r+l) / (r-l),
                            0,           2 / (t-b), 0,          -(t+b) / (t-b),
                            0,           0,         2 / (n-f),  -(n+f) / (n-f),
                            0,           0,         0,          1;

    //TODO: setup projection matrix

    uniform.projective <<   1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 1, 1;

    if (is_perspective)
    {
        // I cannot get perspective to work properly.
        // What i do know about my code:
        //  1) when toggling perspective on, the bunny disappears unless i change 
        //      -e.z() to e.z() in my uniform.cam equation
        //  2) when changing to e.z() my bunny disappears in orthogonal view.
        //  3) when changing to e.z(), my bunny is then backwards in perspective view.
        //      this tells me the camera is facing the wrong way
        //
        // disregard, IDKY but inverting the z-values solved my issue? 
        // any clarity as to whether this was the intended solution or if I have an
        // evident error elsewhere put into the assignment feedback would be appreciated
        
        //TODO setup prespective camera
        // perspective matrix... kinda, Mperspective
        uniform.perspective << n, 0, 0,     0,
                               0, n, 0,     0,
                               0, 0, -(n+f), -(f*n),
                               0, 0, -1,     0;
        //cout << "perspective" << endl;
       
        Matrix4d Mper;
        Mper << uniform.orthographic * uniform.perspective;
        uniform.view =  Mper * uniform.cam;
    }
    else
    {
        uniform.view = uniform.orthographic * uniform.cam;
        //cout << uniform.view << endl << endl;
    }
}

void simple_render(Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) 
    {
        VertexAttributes out;
        out.position = uniform.view * va.position;
        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) 
    {
        
        FragmentAttributes out(uniform.color(0), uniform.color(1), uniform.color(2), uniform.color(3));
        return out;
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) 
    {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    vector<VertexAttributes> vertex_attributes;
    vector<VertexAttributes> thevertices;    
    
    //TODO: build the vertex attributes from vertices and facets
    
    for (int i = 0; i < vertices.size()/3; i++) 
    {
        
        // Create the vertices from the input file
        thevertices.push_back(VertexAttributes(vertices.coeff(i, 0), 
                                            vertices.coeff(i, 1), 
                                            vertices.coeff(i, 2)));
        // cout << vertices.coeff(i, 0) << " " << vertices.coeff(i, 1) << " " << vertices.coeff(i, 2) << endl;
    }
    
    for (int i = 0; i < facets.size()/3; i++) 
    {
        
        // organize the vertices into their triangles
        vertex_attributes.push_back(thevertices[facets.coeff(i, 0)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 1)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 2)]);
        // cout << facets.coeff(i, 0) << " " << facets.coeff(i+facets.size()/3, 0) << " " << facets.coeff(i+ 2 * facets.size()/3, 0) << endl;
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);

}

Matrix4d compute_rotation(const double alpha)
{
    //TODO: Compute the rotation matrix of angle alpha on the y axis around the object barycenter
    // rotation matrix around the y-axis, with angles in radians
    Matrix4d res;
    res <<  cos(alpha * (3.14 / 180)),  0, sin(alpha * (3.14 / 180)),   0,
            0,           1, 0,            0, 
            -sin(alpha * (3.14 / 180)), 0, cos(alpha * (3.14 / 180)),   0,
            0,           0, 0,            1;
    return res;
}

void wireframe_render(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    //compute rotation for GIFs
    Matrix4d trafo = compute_rotation(alpha);
    uniform.view = uniform.view * trafo;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) 
    {
        VertexAttributes out;
        out.position = uniform.view * va.position;
        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) 
    {
        //TODO: fill the shader
        FragmentAttributes out(uniform.color(0), uniform.color(1), uniform.color(2), uniform.color(3));
        return out;
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) 
    {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    vector<VertexAttributes> vertex_attributes;
    vector<VertexAttributes> thevertices;    

    //TODO: generate the vertex attributes for the edges and rasterize the lines

    for (int i = 0; i < vertices.size()/3; i++)
    {
        
        thevertices.push_back(VertexAttributes(vertices.coeff(i, 0), 
                                            vertices.coeff(i, 1), 
                                            vertices.coeff(i, 2)));
    }

    for (int i = 0; i < facets.size()/3; i++) 
    {
        
        vertex_attributes.push_back(thevertices[facets.coeff(i, 0)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 1)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 1)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 2)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 2)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 0)]);
    } 

    //TODO: use the transformation matrix
    // already done above, when rotation is computed

    rasterize_lines(program, uniform, vertex_attributes, 0.5, frameBuffer);
}

void get_shading_program(Program &program) 
{

    program.VertexShader = [](const VertexAttributes &va, 
                                const UniformAttributes &uniform) 
    {

        VertexAttributes out;
        out.position = uniform.view * va.position;
        
        // basically copy and pasted from A3
        Vector3d ray_origin = camera_position;
        Vector3d ray_direction = va.xyz - camera_position;

        Vector3d p = va.xyz;
        Vector3d N = va.normal;

        Vector3d lights_color(0, 0, 0);
        Vector3d color(0, 0, 0);
        for (int i = 0; i < light_positions.size(); ++i)
        {
            const Vector3d light_position = light_positions[i];
            const Vector3d light_color = light_colors[i];

            const Vector3d Li = (light_position - p).normalized();

        Vector3d vector_v = (camera_position - N).normalized();
        Vector3d vector_h = (vector_v + Li).normalized();

        double specular_a = std::pow(std::max(0., N.dot(vector_h)), obj_specular_exponent);
        // Diffuse contribution
        const Vector3d diffuse = obj_diffuse_color * std::max(Li.dot(N), 0.0);

        // Specular contribution, use obj_specular_color
        const Vector3d specular = specular_a * obj_specular_color;

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        
        // light equation
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
        color = ambient_light + lights_color;
        out.color = Vector4d(color.x(), color.y(), color.z(), 1);

        }

        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, 
                                const UniformAttributes &uniform) 
    {

        FragmentAttributes out(va.color[0], va.color[1], va.color[2], va.color[3]);
        out.position = va.position;
        
        return out;    
        };

    program.BlendingShader = [](const FragmentAttributes &fa, 
                            const FrameBufferAttributes &previous) 
    {
        // depth check
        if (fa.position[2] > previous.depth)
        { // if new vertex is closer than previous vertex
            FrameBufferAttributes out(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
            out.depth = fa.position[2];
            return out;
        }
        else
            return previous;    
     };
}

void flat_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, 
                            Eigen::Dynamic, Eigen::Dynamic> &frameBuffer){
    
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    //compute rotation for GIFs
    Eigen::Matrix4d trafo = compute_rotation(alpha);
    uniform.view = uniform.view * trafo;

    vector<VertexAttributes> vertex_attributes;
    vector<VertexAttributes> thevertices; 


    //TODO: generate the vertex attributes for the edges and rasterize the lines

    for (int i = 0; i < vertices.size()/3; i++)
    {
        
        thevertices.push_back(VertexAttributes(vertices.coeff(i, 0), 
                                            vertices.coeff(i, 1), 
                                            vertices.coeff(i, 2)));
    }
    
    for (int i = 0; i < facets.size() / 3; i++) 
    {
        // Get vertices of the triangle
        VertexAttributes v0 = thevertices[facets.coeff(i, 0)];
        VertexAttributes v1 = thevertices[facets.coeff(i, 1)];
        VertexAttributes v2 = thevertices[facets.coeff(i, 2)];    

        // Compute the two vectors on the triangle's plane
        Eigen::Vector3d edge1 = v1.xyz - v0.xyz;
        
        Eigen::Vector3d edge2 = v2.xyz - v0.xyz;

        // Compute the normal vector using cross product
        Vector3d normal = edge1.cross(edge2).normalized();

        // Assign the normal to each vertex of the triangle
        v0.normal = normal;
        v1.normal = normal;
        v2.normal = normal;

        vertex_attributes.push_back(v0);
        vertex_attributes.push_back(v1);
        vertex_attributes.push_back(v2);
    }

    //TODO: set material colors
    //I'll be honest IDK what you mean by this, but output looks right

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

void pv_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    //compute rotation for GIFs
    Eigen::Matrix4d trafo = compute_rotation(alpha);
    uniform.view = uniform.view * trafo;

    vector<VertexAttributes> vertex_attributes;
    vector<VertexAttributes> thevertices; 


    //TODO: generate the vertex attributes for the edges and rasterize the lines

    for (int i = 0; i < vertices.size()/3; i++)
    {
        
        thevertices.push_back(VertexAttributes(vertices.coeff(i, 0), 
                                            vertices.coeff(i, 1), 
                                            vertices.coeff(i, 2)));
    }
    
    for (int i = 0; i < facets.size() / 3; i++) 
    {
        
        // Get vertices of the triangle
        VertexAttributes v0 = thevertices[facets.coeff(i, 0)];
        VertexAttributes v1 = thevertices[facets.coeff(i, 1)];
        VertexAttributes v2 = thevertices[facets.coeff(i, 2)];    

        // Compute the two vectors on the triangle's plane
        Eigen::Vector3d edge1 = v1.xyz - v0.xyz;
        
        Eigen::Vector3d edge2 = v2.xyz - v0.xyz;

        // Compute the normal vector using cross product
        Vector3d normal = edge1.cross(edge2);

        // // Assign the normal to each vertex of the triangle
        thevertices[facets.coeff(i, 0)].normal += normal;
        thevertices[facets.coeff(i, 0)].triangles_touching++;
        thevertices[facets.coeff(i, 1)].normal += normal;
        thevertices[facets.coeff(i, 1)].triangles_touching++;
        thevertices[facets.coeff(i, 2)].normal += normal;
        thevertices[facets.coeff(i, 2)].triangles_touching++;

        
    }
    
    for (int i = 0; i < facets.size() / 3; i++) 
    { // for each triangle
        for (int j = 0; j < 3; j++)
        { // for each vertex in triangle

            if (thevertices[facets.coeff(i, j)].triangles_touching != 0)
            { // check if vertex already normalized, if not, normalize
                thevertices[facets.coeff(i, j)].normal = (thevertices[facets.coeff(i, j)].normal / 
                                    thevertices[facets.coeff(i, j)].triangles_touching).normalized();
                thevertices[facets.coeff(i, j)].triangles_touching = 0; // reset value so when the vertex comes up later the vertex does not get re-normalized
            }
        }
        vertex_attributes.push_back(thevertices[facets.coeff(i, 0)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 1)]);
        vertex_attributes.push_back(thevertices[facets.coeff(i, 2)]);
    }

    //TODO: set material colors
    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

int main(int argc, char *argv[])
{
    setup_scene();

    int W = H * aspect_ratio;
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(W, H);
    vector<uint8_t> image;

    simple_render(frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("simple.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer.setConstant(FrameBufferAttributes()); // reset frame buffer
    
    wireframe_render(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("wireframe.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer.setConstant(FrameBufferAttributes()); // reset frame buffer

    flat_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("flat_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer.setConstant(FrameBufferAttributes()); // reset frame buffer

    pv_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("pv_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);
    frameBuffer.setConstant(FrameBufferAttributes()); // reset frame buffer

    //TODO: add the animation
    
    const char *fileName1 = "wireframe.gif";
    const char *fileName2 = "flat_shading.gif";
    const char *fileName3 = "pv_shading.gif";
    int delay = 25;
    GifWriter w;//wireframe
    GifWriter f; //flat shading
    GifWriter pv; // per-vertex shading

    // Wireframe gif
    GifBegin(&w, fileName1, frameBuffer.rows(), frameBuffer.cols(), delay);
    GifBegin(&f, fileName2, frameBuffer.rows(), frameBuffer.cols(), delay);
    GifBegin(&pv, fileName3, frameBuffer.rows(), frameBuffer.cols(), delay);

    for (double i = 0; i < 360; i += 15)
    {
        cout << "gif creation percentage: " << (i / 360) * 100 << "%" << endl;
        frameBuffer.setConstant(FrameBufferAttributes());
        wireframe_render(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&w, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
        
        frameBuffer.setConstant(FrameBufferAttributes());
        flat_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&f, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
        
        frameBuffer.setConstant(FrameBufferAttributes());
        pv_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&pv, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    cout << "gif creation complete!" << endl;

    GifEnd(&w);
    GifEnd(&f);
    GifEnd(&pv);
    
    return 0;
}