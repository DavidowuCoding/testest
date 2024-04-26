

#include <limits>

#include "Scene.h"

#include "task1.h"

constexpr float epsilon = 0.001f;

bool intersectRayPlane(const float3 &p, const float3 &d, const float4 &plane, float &t)
{
  // TODO implement the intersection test between a ray and a plane.
  //
  // A plane is defined by the plane normal n and the offset w along the normal.
  // The plane is given as parameter plane where (plane.x, plane.y, plane.z) represents the plane normal
  // and plane.w is the offset w.
  //
  // If there is no intersection (Hint: or one we do not care about), return false.
  // Otherwise, compute and set the parameter t such that p + t * d yields the intersection point and return true.

  return false;
}

bool intersectsRayPlane(const float3 &p, const float3 &d, const Plane *planes, std::size_t num_planes, float t_min, float t_max)
{
  // TODO: implement intersection test between a ray and a set of planes.
  // This method only has to detect whether there is an intersection with ANY
  // plane along the given subset of the ray. The ray is given by its start
  // point p and direction d.
  // A plane is defined by the plane normal n and the offset w along the normal.
  // Each plane in planes contains a float4 member p where the plane normal n is
  // (p.x, p.y, p.z) and w is p.w.
  // If an intersection is found that falls on a point on the ray
  // between t_min and t_max, return true. Otherwise, return false.

  return false;
}

const Plane *findClosestHitPlanes(const float3 &p, const float3 &d, const Plane *planes, std::size_t num_planes, float &t)
{
  // TODO: implement intersection test between a ray and a set of planes.
  // This function should find the CLOSEST intersection with a plane along
  // the ray. The ray is given by its start point p and direction d.
  // A plane is defined by the plane normal n and the offset w along the normal.
  // Each plane in planes contains a float4 member p where the plane normal n is
  // (p.x, p.y, p.z) and w is p.w.
  // If an intersection is found, set t to the ray parameter and
  // return a pointer to the hit plane.
  // If no intersection is found, return nullptr.

  return nullptr;
}

bool intersectRayTriangle(const float3 &p, const float3 &d, const float3 &p1,
                          const float3 &p2, const float3 &p3, float &t,
                          float &lambda_1, float &lambda_2)
{
  // TODO implement the intersection test between a ray and a triangle.
  //
  // The triangle is defined by the three vertices p1, p2 and p3
  //
  // If there is no intersection return false.
  // Otherwise, compute and set the parameters lambda_1 and lambda_2
  // to the barycentric coordinates corresponding to the
  // closest point of intersection
  return false;
}
const Triangle *findClosestHitTriangles(const float3 &p, const float3 &d,
                               const Triangle *triangles,
                               std::size_t num_triangles,
                               const float3 *vertices, float &t,
                               float &lambda_1, float &lambda_2)
{
  // TODO: implement intersection test between a ray and a set of triangles.
  // This function should find the CLOSEST intersection with a triangle along
  // the ray. The ray is given by its start point p and direction d. A triangle
  // is represented as an array of three vertex indices. The position of each
  // vertex can be looked up from the vertex array via the vertex index.
  // triangles points to the first element of an array of num_triangles
  // triangles. If an intersection is found, set t to the ray parameter and
  // lambda_1 and lambda_2 to the barycentric coordinates corresponding to the
  // closest point of intersection, and return a pointer to the hit triangle.
  // If no intersection is found, return nullptr.

  return nullptr;
}

bool intersectsRayTriangle(const float3 &p, const float3 &d,
                           const Triangle *triangles, std::size_t num_triangles,
                           const float3 *vertices, float t_min, float t_max)
{
  // TODO: implement intersection test between a ray and a set of triangles.
  // This method only has to detect whether there is an intersection with ANY
  // triangle along the given subset of the ray. The ray is given by its start
  // point p and direction d. A triangle is represented as an array of three
  // vertex indices. The position of each vertex can be looked up from the
  // vertex array via the vertex index. triangles points to an array of
  // num_triangles. If an intersection is found that falls on a point on the ray
  // between t_min and t_max, return true. Otherwise, return false.

  return false;
}

float3 shade(const float3 &p, const float3 &d, const HitPoint &hit,
             const Scene &scene, const Pointlight *lights,
             std::size_t num_lights)
{
  // TODO: implement phong shading.
  // hit represents the surface point to be shaded.
  // hit.position, hit.normal, and hit.k_d and hit.k_s contain the position,
  // surface normal, and diffuse and specular reflection coefficients,
  // hit.m the specular power.
  // lights is a pointer to the first element of an array of num_lights
  // point light sources to consider.
  // Each light contains a member to give its position and color.
  // Return the shaded color.

  // To implement shadows, use scene.intersectsRay(p, d, t_min, t_max) to test
  // whether a ray given by start point p and direction d intersects any
  // object on the section between t_min and t_max.

  return hit.k_d;
}

void render(image2D<float3> &framebuffer, int left, int top, int right,
            int bottom, const Scene &scene, const Camera &camera,
            const Pointlight *lights, std::size_t num_lights,
            const float3 &background_color, int max_bounces)
{
  // TODO: implement raytracing, render the given portion of the framebuffer.
  // left, top, right, and bottom specify the part of the image to compute
  // (inclusive left, top and exclusive right and bottom).
  // camera.eye, camera.lookat, and camera.up specify the position and
  // orientation of the camera, camera.w_s the width of the image plane,
  // and camera.f the focal length to use.
  // Use scene.findClosestHit(p, d) to find the closest point where the ray
  // hits an object.
  // The method returns an std::optional<HitPoint>.
  // If an object is hit, call the function shade to compute the color value
  // of the HitPoint illuminated by the given array of lights.
  // If the ray does not hit an object, use background_color.

  // BONUS: extend your implementation to recursive ray tracing.
  // max_bounces specifies the maximum number of secondary rays to trace.

  for (int y = top; y < bottom; ++y)
  {
    for (int x = left; x < right; ++x)
    {
      framebuffer(x, y) = background_color;
    }
  }
}
