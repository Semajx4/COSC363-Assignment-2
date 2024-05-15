/*==================================================================================
* COSC 363  Computer Graphics
* Department of Computer Science and Software Engineering, University of Canterbury.
*
* A basic ray tracer
* See Lab07.pdf   for details.
*===================================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include <GL/freeglut.h>
#include "Plane.h"
#include "TextureBMP.h"
using namespace std;

const float EDIST = 40.0;
const int NUMDIV = 500;
const int MAX_STEPS = 5;
const float XMIN = -10.0;
const float XMAX = 10.0;
const float YMIN = -10.0;
const float YMAX = 10.0;

const float XMIN_BOX = -50.0;
const float XMAX_BOX = 50.0;
const float YMIN_BOX = -50.0;
const float YMAX_BOX = 50.0;
const float ZMIN_BOX = 0.0;
const float ZMAX_BOX = 400.0;

TextureBMP texture;

vector<SceneObject*> sceneObjects;


//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0);						//Background colour = (0,0,0)
	glm::vec3 lightPos(10, 40, -3);					//Light's position
	glm::vec3 color(0);
	SceneObject* obj;

    ray.closestPt(sceneObjects);					//Compare the ray with all objects in the scene
    if(ray.index == -1) return backgroundCol;		//no intersection
	obj = sceneObjects[ray.index];					//object on which the closest point of intersection is found
	if(ray.index == 0) {
		


		int stripeWidth = 20;
		int ix = (ray.hit.x) / stripeWidth; // Calculate the stripe index in x-direction
		int iz = (ray.hit.z) / stripeWidth; // Calculate the stripe index in z-direction
		int k = (ix + iz) % 2; // Sum of the stripe indices modulo 2 to create a checkered pattern

		if (k == 0) {
			color = glm::vec3(0, 1, 0); // Color 1 for one set of squares
		} else {
			color = glm::vec3(1, 1, 0.5); // Color 2 for the other set of squares
		}
		obj->setColor(color);

	
	}


	color = obj->lighting(lightPos, -ray.dir,ray.hit);						//Object's colour

	glm::vec3 lightVec = lightPos - ray.hit;
	Ray shadowRay(ray.hit, lightVec);
	shadowRay.closestPt(sceneObjects);

	float lightDist = glm::length(lightVec);
	if((shadowRay.index > -1)&&(shadowRay.dist < lightDist)){
		color = 0.2f * obj->getColor();
	}
	
	if(obj->isReflective() && step < MAX_STEPS) {
		float rho = obj->getReflectionCoeff();
		glm::vec3 normalVec = obj->normal(ray.hit);
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
		Ray reflectedRay(ray.hit, reflectedDir);
		glm::vec3 reflectedColor = trace(reflectedRay, step +1);
		color = color + (rho*reflectedColor);
	}

	

	return color;
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX - XMIN) / NUMDIV;  //cell width
	float cellY = (YMAX - YMIN) / NUMDIV;  //cell height
	glm::vec3 eye(0., 0., 0.);

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a tiny quad.

	for (int i = 0; i < NUMDIV; i++)	//Scan every cell of the image plane
	{
		xp = XMIN + i * cellX;
		for (int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j * cellY;

			glm::vec3 dir(xp + 0.5 * cellX, yp + 0.5 * cellY, EDIST);	//direction of the primary ray

			Ray ray = Ray(eye, dir);

			glm::vec3 col = trace(ray, 1); //Trace the primary ray and get the colour value
			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp + cellX, yp);
			glVertex2f(xp + cellX, yp + cellY);
			glVertex2f(xp, yp + cellY);
		}
	}

    glEnd();
    glFlush();
}



//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL 2D orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);

    glClearColor(0, 0, 0, 1);

	texture = TextureBMP("/home/james/Documents/Assignment-2/363-assignment-2/Butterfly.bmp");

	// Define the planes of the Cornell box
	Plane *floorPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX),    // Point A
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX),   // Point B
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX),    // Point C
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX)   // Point D
		
		
	);
	floorPlane->setColor(glm::vec3(0.5, 0.8, 0));   //Set colour 
	sceneObjects.push_back(floorPlane);

	Plane *ceilingPlane = new Plane(
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX),    // Point A
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX),    // Point B
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX),    // Point C
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX)     // Point D
	);
	ceilingPlane->setColor(glm::vec3(0.58, 0.85, 0));   //Set colour 
	sceneObjects.push_back(ceilingPlane);

	Plane *leftWallPlane = new Plane(
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX),    // Point A
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX),    // Point B
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX),    // Point C
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX)     // Point D
	);
	leftWallPlane->setColor(glm::vec3(0.9, 0.8, 0));   //Set colour
	sceneObjects.push_back(leftWallPlane);

	Plane *rightWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX),    // Point A
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX),    // Point B
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX),    // Point C
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX)     // Point D
	);
	rightWallPlane->setColor(glm::vec3(0.5, 0.4, 0));   //Set colour
	sceneObjects.push_back(rightWallPlane);

	Plane *backWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX),    // Point A
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX),    // Point B
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX),    // Point C
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX)     // Point D
	);
	backWallPlane->setColor(glm::vec3(0.2, 0.8, 0));   //Set colour
	backWallPlane->setReflectivity(false);
	sceneObjects.push_back(backWallPlane);




	Sphere *sphere1 = new Sphere(glm::vec3(-5.0, 0.0, 150.0), 15.0);
	sphere1->setColor(glm::vec3(0, 0, 1));   //Set colour to blue
	sphere1->setReflectivity(true, 0.8);
	sceneObjects.push_back(sphere1);		 //Add sphere to scene objects

	Sphere *sphere2 = new Sphere(glm::vec3(5.0, -10.0, 140.0), 5.0);
	sphere2->setColor(glm::vec3(0, 1, 1));   //Set colour to blue
	sceneObjects.push_back(sphere2);		 //Add sphere to scene objects

	Sphere *sphere3 = new Sphere(glm::vec3(5.0, 5.0, 90.0), 4.0);
	sphere3->setColor(glm::vec3(1, 0, 1));   //Set colour to blue
	sceneObjects.push_back(sphere3);		 //Add sphere to scene objects

	Sphere *sphere4 = new Sphere(glm::vec3(10.0, 10.0, 105.0), 3.0);
	sphere4->setColor(glm::vec3(1, 0, 0));   //Set colour to blue
	sceneObjects.push_back(sphere4);		 //Add sphere to scene objects



}


int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracing");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
