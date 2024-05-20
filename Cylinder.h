#ifndef H_CYLINDER
#define H_CYLINDER
#include <glm/glm.hpp>
#include "SceneObject.h"

/**
 * Defines a simple cylinder located at 'center'
 * with specified radius and specified height
*/
class Cylinder : public SceneObject
{
private:
    glm::vec3 center = glm::vec3(0);
    float radius = 1;
    float height = 1;
public:
    Cylinder() {};

    Cylinder(glm::vec3 c, float r, float h) : center(c), radius(r), height(h) {};
    
    float intersect(glm::vec3 p0, glm::vec3 dir);
    
    glm::vec3 normal(glm::vec3 point);

};

#endif //!H_SPHERE
