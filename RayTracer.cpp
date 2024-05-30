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
#include "Cylinder.h"
#include "SceneObject.h"
#include "Ray.h"
#include <GL/freeglut.h>
#include "Plane.h"
#include "TextureBMP.h"
#include <random>
#include <glm/gtc/matrix_transform.hpp>
using namespace std;
const bool FOGGY = false;
const bool ANTI_ALIAS = true;
const float EDIST = 50.0;
const int NUMDIV = 1000;
const int MAX_STEPS = 5;
const int SS_FACTOR = 2;
const float XMIN = -20.0;
const float XMAX = 20.0;
const float YMIN = -25.0;
const float YMAX = 25.0;

const float XMIN_BOX = -50.0;
const float XMAX_BOX = 50.0;
const float YMIN_BOX = -50.0;
const float YMAX_BOX = 50.0;
const float ZMIN_BOX = -50.0;
const float ZMAX_BOX = 50.0;
TextureBMP texture;

TextureBMP earthTex;

vector<SceneObject *> sceneObjects;

float getShadowFactor(SceneObject *obj)
{
	if (obj->isTransparent())
	{
		return 0.4f; // Lighter shadow for transparent objects
	}
	if (obj->isRefractive())
	{
		return 0.4f; // Lighter shadow for transparent objects
	}
	else
	{
		return 0.2f; // Darker shadow for opaque objects
	}
}


//---The most important function in a ray tracer! ----------------------------------
//   Computes the colour value obtained by tracing a ray and finding its
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0);			// Background colour = (0,0,0)
	glm::vec3 lightPos1(0., 45., -30.); // Light's position
	glm::vec3 color(0.5, 0.5, 0.2);
	SceneObject *obj;

	ray.closestPt(sceneObjects); // Compare the ray with all objects in the scene
	if (ray.index == -1)
		return backgroundCol;	   // no intersection
	obj = sceneObjects[ray.index]; // object on which the closest point of intersection is found
	if (ray.index == 0)
	{
		float squareSize = 50;
		bool isZ = glm::mod(ray.hit.z, squareSize) > squareSize * 0.5;
		bool isX = glm::mod(ray.hit.x, squareSize) > squareSize * 0.5;
		if (isX xor isZ)
		{
			color = glm::vec3(0.5, 0.5, 0.5); //Color grey for first set of squares
		}
		else
		{
			color = glm::vec3(1, 1, 1); // Color white for the other set of squares
		}
		obj->setColor(color);
	}
	if(ray.index == 1) {
		float blueSize = (30);
		float redSize = (50);

		bool isZBlue = glm::mod(ray.hit.z,blueSize) > blueSize*0.5;
		bool isXBlue = glm::mod(ray.hit.x,blueSize) >blueSize*0.5;
		bool isZRed = glm::mod(ray.hit.z,redSize) >redSize*0.5;
		bool isXRed = glm::mod(ray.hit.x,redSize) >redSize*0.5;

		if ((isXBlue xor isZBlue)){
			color = glm::vec3(0,0,1);
		} else if ((isXRed xor isZRed )) {
			color = glm::vec3(1,0,0);
		} else {
			color = glm::vec3(1.,1.,1.);
		}
		obj->setColor(color);
	} 

	if (ray.index == 5)
	{
		obj->setColor(glm::vec3(1,1, 1));
		float x1 = 40;
		float x2 = 0;
		float y1 = 0;
		float y2 = -40;
		float texcoordt = (ray.hit.x - x1) / (x2 - x1);
		float texcoords = (ray.hit.y - y1) / (y2 - y1);
		if (
			texcoords > 0 &&
			texcoords < 1 &&
			texcoordt > 0 &&
			texcoordt < 1)
		{
			color = texture.getColorAt(texcoords, texcoordt);
			obj->setColor(color);
		}
	}


	if (ray.index == 10 )
	{
		glm::vec3 normalVec = obj->normal(ray.hit);

        float texcoordt = atan2(-normalVec.z, normalVec.x)/(2*M_PI) + 0.5;
        float texcoords =  -asin(normalVec.y)/ M_PI + 0.5;
        if(texcoords > 0 && texcoords < 1 && texcoordt > 0 && texcoordt < 1) {
            color=earthTex.getColorAt(texcoords, texcoordt);
            obj->setColor(color);
        }
	}

	color = obj->lighting(lightPos1, -ray.dir, ray.hit); // Object's colour

	glm::vec3 lightVec = lightPos1 - ray.hit;
	Ray shadowRay(ray.hit, lightVec);
	shadowRay.closestPt(sceneObjects);

	float lightDist = glm::length(lightVec);
	float shadowFactor = 1.0f; // Default shadow factor (no shadow

	if ((shadowRay.index > -1) && (shadowRay.dist < lightDist))
	{

		SceneObject *shadowObj = sceneObjects[shadowRay.index];
		shadowFactor = getShadowFactor(shadowObj);
	}
	color *= shadowFactor;



	if (obj->isReflective() && step < MAX_STEPS)
	{
		float rho = obj->getReflectionCoeff();
		glm::vec3 normalVec = obj->normal(ray.hit);
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
		Ray reflectedRay(ray.hit, reflectedDir);
		glm::vec3 reflectedColor = trace(reflectedRay, step + 1);
		color = color + (rho * reflectedColor);
	}

	if (obj->isTransparent() && step < MAX_STEPS)
	{
		float transparencyCoeff = obj->getTransparencyCoeff(); // Get transparency coefficient
		glm::vec3 normalVec = obj->normal(ray.hit);			   // Normal at the hit point
		glm::vec3 transmittedDir = ray.dir;
		glm::vec3 offset = normalVec * 0.001f;
		Ray transmittedRay(ray.hit + offset, transmittedDir); // Create the transmitted ray
		glm::vec3 transmittedColor = trace(transmittedRay, step + 1);
		color = (1.0f - transparencyCoeff) * color + transparencyCoeff * transmittedColor;
	}

	if (obj->isRefractive() && step < MAX_STEPS)
	{
		float refractiveIndex = obj->getRefractiveIndex();
		glm::vec3 normalVec = obj->normal(ray.hit);
		glm::vec3 transmittedDir = ray.dir;
		glm::vec3 refractedDirection;

		float cosineOfOIncidentAngle = glm::dot(-transmittedDir, normalVec);
		float refractiveIndexOutside = 1.0f, refractiveIndexInside = refractiveIndex;

		if (cosineOfOIncidentAngle < 0)
		{
			cosineOfOIncidentAngle = -cosineOfOIncidentAngle;
			std::swap(refractiveIndexOutside, refractiveIndexInside);
			normalVec = -normalVec;
		}
		float etaRatio = refractiveIndexOutside / refractiveIndexInside;
		float refractionFactor = 1 - etaRatio * etaRatio * (1 - cosineOfOIncidentAngle * cosineOfOIncidentAngle);
		if (refractionFactor < 0)
		{
			refractedDirection = glm::reflect(ray.dir, normalVec);
		}
		else
		{
			refractedDirection = etaRatio * ray.dir + (etaRatio * cosineOfOIncidentAngle - sqrtf(refractionFactor)) * normalVec;
		}
		glm::vec3 offset = normalVec * 0.001f;
		Ray refractedRay(ray.hit + offset, refractedDirection);

		
		glm::vec3 refractedColor = trace(refractedRay, step + 1);

		
		float refractionCoeff = obj->getRefractionCoeff();
		color = (1.0f - refractionCoeff) * color + refractionCoeff * refractedColor;
	}
	if (FOGGY){
		float z1 = ZMIN_BOX;
		float z2 = ZMAX_BOX+50; //most foggy should be past the back of the box so that it is not completly white there.
		float lambda = ((ray.hit.z)-z1)/(z2-z1);
		return (1-lambda)*color + lambda*glm::vec3(1,1,1);
	} else {
		return color;
	}
	
}


/*
Function to determine if a pixel is the edge of two or more color values
*/
bool isEdge(int i, int j, const std::vector<std::vector<glm::vec3>>& colors, float threshold = 0.01f) {
    glm::vec3 targetColor = colors[i][j];
    std::vector<std::pair<int, int>> offsets = {
        {-1, -1}, {-1, 0}, {-1, 1},
        {0, -1},           {0, 1},
        {1, -1},  {1, 0},  {1, 1}
    };

    for (const auto& offset : offsets) {
        int ni = i + offset.first;
        int nj = j + offset.second;
        if (ni >= 0 && ni < colors.size() && nj >= 0 && nj < colors[0].size()) {
            glm::vec3 neighborColor = colors[ni][nj];
            float colorDifference = glm::length(targetColor - neighborColor);
            if (colorDifference > threshold) {
                return true;
            }
        }
    }
    return false;
}


//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;						  // grid point
	float cellX = (XMAX - XMIN) / NUMDIV; // cell width
	float cellY = (YMAX - YMIN) / NUMDIV; // cell height
	glm::vec3 eye(0., 0., -200.);

	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_QUADS); // Each cell is a tiny quad.

	std::vector<std::vector<glm::vec3>> colors(NUMDIV, std::vector<glm::vec3>(NUMDIV));

	for (int i = 0; i < NUMDIV; i++) // Scan every cell of the image plane in order to store the color values of each cell
	{
		xp = XMIN + i * cellX;
		for (int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j * cellY;
			glm::vec3 dir(xp + 0.5 * cellX, yp + 0.5 * cellY, EDIST);	//direction of the primary ray
			Ray ray = Ray(eye, dir);
			glm::vec3 col = trace(ray, 1); //Trace the primary ray and get the colour value
			colors[i][j] = col; //set colors i,j to its color val.
		}
	}

	for (int i = 0; i < NUMDIV; i++) { //Loop through all the cells again
        xp = XMIN + i * cellX;
        for (int j = 0; j < NUMDIV; j++) {
            yp = YMIN + j * cellY;

            glm::vec3 cellColour = colors[i][j];
            if (isEdge(i, j, colors) && ANTI_ALIAS) { //Check if the cell is and edge
                glm::vec3 col_avg(0.0f);
                for (int k = 0; k < SS_FACTOR; k++) {
                    for (int l = 0; l < SS_FACTOR; l++) {
                        glm::vec3 dir(xp + k * cellX / SS_FACTOR + 0.5 * cellX / SS_FACTOR, yp + l * cellY / SS_FACTOR + 0.5 * cellY / SS_FACTOR, EDIST);
                        Ray ray = Ray(eye, dir);
                        glm::vec3 col = trace(ray, 1);
                        col_avg += col;
                    }
                }
                col_avg /= (SS_FACTOR * SS_FACTOR);
                cellColour = col_avg;
            }

            glColor3f(cellColour.r, cellColour.g, cellColour.b);
            glVertex2f(xp, yp);
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

	texture = TextureBMP("Butterfly.bmp"); //Butterfly texture
	earthTex = TextureBMP("Earth.bmp"); //Earth texture

	// Define the planes of the Cornell box
	Plane *floorPlane = new Plane( //Index 0 -> Will be colored in a checkered pattern
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX), 
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX), 
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX), 
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX)	

	);
	sceneObjects.push_back(floorPlane);

	Plane *ceilingPlane = new Plane( //Index 1 -> will be colored in a procedurally generated pattern
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX), 
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX), 
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX), 
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX)	 
	);
	sceneObjects.push_back(ceilingPlane);

	Plane *rightWallPlane = new Plane(
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX), 
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX), 
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX), 
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX)	 
	);
	rightWallPlane->setColor(glm::vec3(0.2, 0.8, 0.2)); // Set colour: Green
	sceneObjects.push_back(rightWallPlane);

	Plane *leftWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX), 
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX), 
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX), 
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX)	 
	);
	leftWallPlane->setColor(glm::vec3(0.7, 0.2, 0.2)); // Set colour: Red
	sceneObjects.push_back(leftWallPlane);

	Plane *backWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX), 
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX), 
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX), 
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX)	 
	);
	backWallPlane->setColor(glm::vec3(0.0, 0.0, 0.0)); // Set colour: black
	backWallPlane->setReflectivity(true, 1); //Make surface completly reflective.
	sceneObjects.push_back(backWallPlane);

	Plane *frontWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX), // Point B
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX), // Point C
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX), // Point D
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX)	 // Point A
	);
	frontWallPlane->setColor(glm::vec3(0, 0, 1)); // Set colour
	sceneObjects.push_back(frontWallPlane);

	glm::mat4 transform = glm::mat4(1.0); //Set transformation matrix to the identity matrix
	transform = glm::scale(transform, glm::vec3(1, 1.1, 1)); //Scale in y direction by 1.1

	Sphere *transformedSphere = new Sphere(glm::vec3(-25.0, -20.0, 20.0), 10.0);
	transformedSphere->setColor(glm::vec3(0, 0, 0)); // Set color: black
	transformedSphere->setReflectivity(true, .8); //Set to reflective
	transformedSphere->setTransform(true, transform); //Squish the sphere
	sceneObjects.push_back(transformedSphere);

	Sphere *refractiveSphere = new Sphere(glm::vec3(25.0, -35, -40.0), 8.0);
	refractiveSphere->setColor(glm::vec3(1, 1, 1)); //set color: white
	refractiveSphere->setRefractivity(true, 1, 1.5);
	sceneObjects.push_back(refractiveSphere);

	Sphere *transparentSphere = new Sphere(glm::vec3(-20.0, -35, -30.0), 8.0);
	transparentSphere->setColor(glm::vec3(1, 1, 1)); // set color: black
	transparentSphere->setTransparency(true, 0.8);
	transparentSphere->setReflectivity(true, 0.1);
	sceneObjects.push_back(transparentSphere);

	Cylinder *cyl = new Cylinder(glm::vec3(15, -50, 25.0),2, 20);
	cyl->setColor(glm::vec3(0, 0, 1)); // Set colour to blue
	sceneObjects.push_back(cyl);
	
	Sphere *texSphere = new Sphere(glm::vec3(35.0, -20, 30.0), 12.0);
	texSphere->setColor(glm::vec3(1, 1, 1));  
	sceneObjects.push_back(texSphere);
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(1000, 1000);
	glutInitWindowPosition(20, 20);
	glutCreateWindow("Raytracing");

	glutDisplayFunc(display);
	initialize();

	glutMainLoop();
	return 0;
}
