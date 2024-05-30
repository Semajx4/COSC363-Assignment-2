/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The sphere class
*  This is a subclass of SceneObject, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Sphere.h"
#include <math.h>

/**
* Sphere's intersection method.  The input is a ray. 
*/
float Sphere::intersect(glm::vec3 p0, glm::vec3 dir)
{
    glm::mat4 transform = getTransform();
    p0 = transform*glm::vec4(p0,1.0);
    dir = transform*glm::vec4(dir,1.0);
    float tScale = glm::length(dir);
    dir = glm::normalize(dir);
    glm::vec3 vdif = p0 - center;   //Vector s (see Slide 28)
    float b = glm::dot(dir, vdif);
    float len = glm::length(vdif);
    float c = len*len - radius*radius;
    float delta = b*b - c;
   
	if(delta < 0.001) return -1.0;    //includes zero and negative values

    float t1 = -b - sqrt(delta);
    float t2 = -b + sqrt(delta);
    t1 /=tScale;
    t2 /= tScale;

	if (t1 < 0)
	{
		return (t2 > 0) ? t2 : -1;
	}
	else return t1;
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the sphere.
*/
glm::vec3 Sphere::normal(glm::vec3 p)
{
    
    glm::mat4 normalTransform = glm::transpose(glm::inverse(transmat_)); // Inverse transpose of the model matrix

    p = transmat_ * glm::vec4(p, 1.0f);
    glm::vec3 n = p - center;

    // Transform the normal vector
    n = glm::vec3(normalTransform * glm::vec4(n, 0.0f)); // 0.0f as w-component for direction

    n = glm::normalize(n);
    return n;
}
