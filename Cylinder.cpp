/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of SceneObject, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include <math.h>

/**
* Cylinder's intersection method.  The input is a ray. 
*/
float Cylinder::intersect(glm::vec3 p0, glm::vec3 dir){

     // Define the cylinder's axis direction (assuming the cylinder's axis is aligned with the y-axis)
    glm::vec3 cylAxis = glm::vec3(0.0f, 1.0f, 0.0f);

    // Compute the vector from the ray origin to the cylinder center
    glm::vec3 U = p0 - center;
    
    // Components of the direction vector orthogonal and parallel to the cylinder's axis
    glm::vec3 V_perp = dir - glm::dot(dir, cylAxis) * cylAxis;
    glm::vec3 U_perp = U - glm::dot(U, cylAxis) * cylAxis;

    float a = glm::dot(V_perp, V_perp);
    float b = 2.0f * glm::dot(V_perp, U_perp);
    float c = glm::dot(U_perp, U_perp) - radius * radius;

    // Calculate the discriminant
    float discriminant = b * b - 4 * a * c;

    if (discriminant >= 1e-4) {
        float sqrtDiscriminant = sqrt(discriminant);
        float t1 = (-b - sqrtDiscriminant) / (2.0f * a);
        float t2 = (-b + sqrtDiscriminant) / (2.0f * a);

        // Calculate intersection points
        glm::vec3 intersection1 = p0 + t1 * dir;
        glm::vec3 intersection2 = p0 + t2 * dir;

        // Check if the intersection points are within the height bounds of the cylinder
        bool isInBounds1 = (intersection1.y >= center.y) && (intersection1.y <= center.y + height);
        bool isInBounds2 = (intersection2.y >= center.y) && (intersection2.y <= center.y + height);

        // Determine the valid intersection point
        if (t1 > 0.0f && isInBounds1) {
            return t1;
        } else if (t2 > 0.0f && isInBounds2) {
            return t2;
        }

    }
       // Check intersection with the bottom cap
    float tBottom = (center.y - p0.y) / dir.y;
    glm::vec3 bottomIntersection = p0 + tBottom * dir;
    if (tBottom > 0 && glm::distance(glm::vec2(bottomIntersection.x, bottomIntersection.z), glm::vec2(center.x, center.z)) <= radius) {
        return tBottom;
    }

    float tTop = (center.y + height - p0.y) / dir.y;
    glm::vec3 topIntersection = p0 + tTop * dir;
    if (tTop > 0 && glm::distance(glm::vec2(topIntersection.x, topIntersection.z), glm::vec2(center.x, center.z)) <= radius) {
        return tTop;
    }
    // If both intersections are invalid or behind the ray origin, return -1
    return -1.0f;


};


/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the cylinder.
*/
glm::vec3 Cylinder::normal(glm::vec3 p)
{   
    // Check if the point is on the bottom cap
    if (fabs(p.y - center.y) < 1e-4) {
        return glm::vec3(0.0f, -1.0f, 0.0f);
    }
    // Check if the point is on the top cap
    if (fabs(p.y - (center.y + height)) < 1e-4) {
        return glm::vec3(0.0f, 1.0f, 0.0f);
    }
    glm::vec3 n = glm::vec3(p.x - center.x, 0.0f, p.z - center.z);
    n = glm::normalize(n);
    return n;
};