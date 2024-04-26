

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
    float3 n = {plane.x, plane.y, plane.z}; // Extract normal vector from plane
    float w = plane.w; // Distance along the normal from the origin

    float denominator = dot(n, d);

    if (fabs(denominator) > 1e-6) { // Check for non-zero denominator to avoid division by zero
        float numerator = dot(n, p) + w;
        t = -numerator / denominator;
        if (t >= 0) {
            return true; // Intersection occurs if t is non-negative
        }
    }
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
    const Plane *closest_plane = nullptr;
    float closest_t = std::numeric_limits<float>::infinity(); // Initialize with max float to find the minimum
    float current_t;

    for (std::size_t i = 0; i < num_planes; ++i) {
        if (intersectRayPlane(p, d, planes[i].p, current_t)) {
            if (current_t < closest_t) {
                closest_t = current_t;
                closest_plane = &planes[i];
            }
        }
    }

    if (closest_plane) {
        t = closest_t; // Update t to the closest intersection distance
        return closest_plane;
    }

    return nullptr; // Return nullptr if no intersection is found
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
    float3 e1 = p2 - p1; // Edge 1
    float3 e2 = p3 - p1; // Edge 2
    float3 h = cross(d, e2);
    float a = dot(e1, h);

    if (fabs(a) < 1e-6) // This means the ray is parallel to the triangle.
        return false;

    float f = 1.0f / a;
    float3 s = p - p1;
    float u = f * dot(s, h);

    if (u < 0.0f || u > 1.0f)
        return false; // The intersection is outside of the triangle.

    float3 q = cross(s, e1);
    float v = f * dot(d, q);

    if (v < 0.0f || u + v > 1.0f)
        return false; // The intersection is outside of the triangle.

    // Compute the value of t, the scalar distance along the ray to the intersection point.
    t = f * dot(e2, q);
    if (t > epsilon) { // Check against epsilon to avoid intersections very close to the ray origin
        lambda_1 = u;
        lambda_2 = v;
        return true;
    }

    return false; // Either behind the ray's start or too close to it

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
    const Triangle *closest_triangle = nullptr;
    float closest_t = std::numeric_limits<float>::infinity();
    float current_t, current_lambda_1, current_lambda_2;

    for (std::size_t i = 0; i < num_triangles; ++i) {
        const float3& p1 = vertices[triangles[i][0]];
        const float3& p2 = vertices[triangles[i][1]];
        const float3& p3 = vertices[triangles[i][2]];

        if (intersectRayTriangle(p, d, p1, p2, p3, current_t, current_lambda_1, current_lambda_2)) {
            if (current_t < closest_t && current_t > epsilon) { // Ensure we're not too close to the ray origin
                closest_t = current_t;
                lambda_1 = current_lambda_1;
                lambda_2 = current_lambda_2;
                closest_triangle = &triangles[i];
            }
        }
    }

    if (closest_triangle) {
        t = closest_t; // Update t to the closest intersection distance
        return closest_triangle;
    }

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
    // Calculate the camera basis vectors
    float3 w = normalize(camera.lookat - camera.eye);  // World forward vector (z-axis)
    float3 u = normalize(cross(camera.up, w));         // Camera's right vector (x-axis)
    float3 v = cross(w, u);                            // Camera's up vector (y-axis)

    // Calculate image plane dimensions based on camera parameters
    float image_plane_width = camera.w_s;
    float image_plane_height = image_plane_width * (float)(bottom - top) / (right - left);

    // Calculate pixel dimensions on the image plane
    float pixel_width = image_plane_width / (right - left);
    float pixel_height = image_plane_height / (bottom - top);

    for (int y = top; y < bottom; ++y) {
        for (int x = left; x < right; ++x) {
            // Calculate normalized device coordinates (NDC) for the current pixel
            float ndc_x = (x - left + 0.5f) / (right - left) - 0.5f;
            float ndc_y = (y - top + 0.5f) / (bottom - top) - 0.5f;

            // Convert NDC to world coordinates on the image plane
            float3 pixel_world = camera.eye + w * camera.f + u * ndc_x * image_plane_width + v * ndc_y * image_plane_height;

            // Compute the ray direction
            float3 ray_direction = normalize(pixel_world - camera.eye);

            // Ray origin is the camera's position
            float3 ray_origin = camera.eye;

            // Find the closest hit point
            auto hit_point = scene.findClosestHit(ray_origin, ray_direction);

            if (hit_point) {
                // Calculate color via shading
                //framebuffer(x, y) = shade(ray_origin, ray_direction, *hit_point, scene, lights, num_lights);
                framebuffer(x, y) = {0,0,1};
            } else {
                // Set background color if no hit
                framebuffer(x, y) = background_color;
            }
        }
    }
}
