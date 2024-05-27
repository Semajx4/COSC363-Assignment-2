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

const float EDIST = 40.0;
const int NUMDIV = 1000;
const int MAX_STEPS = 5;
const int SS_FACTOR = 1;
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

bool isBetweenWaves(float rayhitx, float rayhity, float stepSize, float amp, float omega, float yshift){
	float sinX = amp*sin((1/omega)*rayhitx)+yshift;
	float sinXwStep = amp*sin((1/omega)*rayhitx )+stepSize+yshift;
	float y = glm::mod(rayhity, (stepSize + omega));
	bool out = (y > sinX) && (y < sinXwStep);
	return out;
}

//---The most important function in a ray tracer! ----------------------------------
//   Computes the colour value obtained by tracing a ray and finding its
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0);			// Background colour = (0,0,0)
	glm::vec3 lightPos1(0., 45., -30.); // Light's position
	// glm::vec3 lightPos2(-30., 100., -15.);					//Light's position
	glm::vec3 color(0.5, 0.5, 0.2);
	SceneObject *obj;

	ray.closestPt(sceneObjects); // Compare the ray with all objects in the scene
	if (ray.index == -1)
		return backgroundCol;	   // no intersection
	obj = sceneObjects[ray.index]; // object on which the closest point of intersection is found
	if (ray.index == 0)
	{
		float squareSize = 10;
		bool isZ = glm::mod(ray.hit.z, squareSize) > squareSize * 0.5;
		bool isX = glm::mod(ray.hit.x, squareSize) > squareSize * 0.5;
		if (isX xor isZ)
		{
			color = glm::vec3(0, 0, 0);
		}
		else
		{
			color = glm::vec3(1, 1, 1); // Color 2 for the other set of squares
		}
		obj->setColor(color);
	}
	if (ray.index == 1)
	{
		float stepSize = 5;
		float amp = 2.50;
		float omega = 5.00;
		float yshift = 2.50;
		
		if (ray.hit.x > 0){
			if (isBetweenWaves(ray.hit.x, ray.hit.y, stepSize,amp, omega, yshift)){
				color = glm::vec3(1,1,1);
			} else {
				color = glm::vec3(0,0,0);
			}

		} else{
			if (isBetweenWaves(-ray.hit.x, ray.hit.y, stepSize,amp, omega, yshift)){
				color = glm::vec3(1,1,1);
			} else {
				color = glm::vec3(0,0,0);
			}


		}
		obj->setColor(color);
	}
	if (ray.index == 5)
	{
		obj->setColor(glm::vec3(1, 1, 1));
		float x1 = -20;
		float x2 = 20;
		float y1 = -20;
		float y2 = 20;
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

	color = obj->lighting(lightPos1, -ray.dir, ray.hit); // Object's colour

	glm::vec3 lightVec = lightPos1 - ray.hit;
	Ray shadowRay(ray.hit, lightVec);
	shadowRay.closestPt(sceneObjects);

	float lightDist = glm::length(lightVec);
	float shadowFactor = 1.0f; // Default shadow factor (no shadow)

	if ((shadowRay.index > -1) && (shadowRay.dist < lightDist))
	{

		SceneObject *shadowObj = sceneObjects[shadowRay.index];
		shadowFactor = getShadowFactor(shadowObj);
	}
	color *= shadowFactor;
	// Lighting calculations for the second light source
	// glm::vec3 colorFromSecondLight = obj->lighting(lightPos2, -ray.dir, ray.hit);
	// glm::vec3 lightVec2 = lightPos2 - ray.hit;
	// Ray shadowRay2(ray.hit, lightVec2);
	// shadowRay2.closestPt(sceneObjects);
	// float lightDist2 = glm::length(lightVec2);
	// if((shadowRay2.index > -1) && (shadowRay2.dist < lightDist2)){
	// colorFromSecondLight = 0.2f * obj->getColor(); // Apply shadow color
	//  }

	// color += colorFromSecondLight;

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
			// Total internal reflection, treat as a perfect mirror reflection
			refractedDirection = glm::reflect(ray.dir, normalVec);
		}
		else
		{
			refractedDirection = etaRatio * ray.dir + (etaRatio * cosineOfOIncidentAngle - sqrtf(refractionFactor)) * normalVec;
		}
		glm::vec3 offset = normalVec * 0.001f;
		Ray refractedRay(ray.hit + offset, refractedDirection); // Create the refracted ray

		// Trace the refracted ray and get its color
		glm::vec3 refractedColor = trace(refractedRay, step + 1);

		// Combine the original color with the refracted color
		float refractionCoeff = obj->getRefractionCoeff();
		color = (1.0f - refractionCoeff) * color + refractionCoeff * refractedColor;
	}

	return color;
}
#include <vector>
#include <glm/glm.hpp>

bool similarColoursNearby(int i, int j, const std::vector<std::vector<glm::vec3>>& colors) {
    // Define the threshold for color similarity
    float threshold = 0.1f; // Adjust as needed
    
    // Get the color of the cell at position (i, j)
    glm::vec3 targetColor = colors[i][j];
    
    // Define the neighboring offsets (assuming 8-connected neighbors)
    std::vector<std::pair<int, int>> offsets = {
        {-1, -1}, {-1, 0}, {-1, 1},
        {0, -1},           {0, 1},
        {1, -1},  {1, 0},  {1, 1}
    };
    
    // Iterate over the neighboring cells
    for (const auto& offset : offsets) {
        int ni = i + offset.first;
        int nj = j + offset.second;
        
        // Check if the neighboring cell is within bounds
        if (ni >= 0 && ni < colors.size() && nj >= 0 && nj < colors[0].size()) {
            // Get the color of the neighboring cell
            glm::vec3 neighborColor = colors[ni][nj];
            
            // Calculate the color difference between the target color and the neighbor color
            float colorDifference = glm::length(targetColor - neighborColor);
            
            // Check if the color difference is within the threshold
            if (colorDifference > threshold) {
                // Colors are not similar, return false
                return false;
            }
        }
    }
    
    // All neighboring colors are similar, return true
    return true;

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

	for (int i = 0; i < NUMDIV; i++) // Scan every cell of the image plane
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

	for (int i = 0; i < NUMDIV; i++) // Scan every cell of the image plane
	{
		xp = XMIN + i * cellX;
		for (int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j * cellY;


			//if adjacent cells are not the same colour
			glm::vec3 cellColour = colors[i][j];
			if(!similarColoursNearby(i,j, colors)){
				glm::vec3 col_avg(0.0f);
				for (int k = 0; k < SS_FACTOR; k++)
				{
					for (int l = 0; l < SS_FACTOR; l++)
					{
						glm::vec3 dir(xp + (k + (float)rand() / RAND_MAX) * cellX,
									yp + (l + (float)rand() / RAND_MAX) * cellY,
									EDIST);
						Ray ray = Ray(eye, dir);
						glm::vec3 col = trace(ray, 1); // Trace the primary ray and get the colour value
						col_avg += col;
					}
				}

				col_avg /= (SS_FACTOR * SS_FACTOR);
				cellColour = col_avg;
			} 

			
			glColor3f(cellColour.r, cellColour.g, cellColour.b);
			glVertex2f(xp, yp); // Draw each cell with its color value
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

	texture = TextureBMP("/csse/users/jsu103/Documents/cosc363/363-assignment-2/363-assignment-2/Butterfly.bmp");

	// Define the planes of the Cornell box
	Plane *floorPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX), // Point A
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX), // Point B
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX), // Point C
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX)	 // Point D

	);
	floorPlane->setColor(glm::vec3(0.5, 0.8, 0)); // Set colour
	sceneObjects.push_back(floorPlane);

	Plane *ceilingPlane = new Plane(
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX), // Point A
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX), // Point B
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX), // Point C
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX)	 // Point D
	);
	ceilingPlane->setColor(glm::vec3(1, 1, 1)); // Set colour
	ceilingPlane->setReflectivity(true, 0.01);
	sceneObjects.push_back(ceilingPlane);

	Plane *rightWallPlane = new Plane(
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX), // Point A
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX), // Point B
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX), // Point C
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX)	 // Point D
	);
	rightWallPlane->setColor(glm::vec3(0, 0.8, 0)); // Set colour
	sceneObjects.push_back(rightWallPlane);

	Plane *leftWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX), // Point A
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX), // Point B
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX), // Point C
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX)	 // Point D
	);
	leftWallPlane->setColor(glm::vec3(0.8, 0, 0)); // Set colour
	sceneObjects.push_back(leftWallPlane);

	Plane *backWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMAX_BOX), // Point A
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMAX_BOX), // Point B
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMAX_BOX), // Point C
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMAX_BOX)	 // Point D
	);
	backWallPlane->setColor(glm::vec3(1, 1, 1)); // Set colour
	sceneObjects.push_back(backWallPlane);

	Plane *frontWallPlane = new Plane(
		glm::vec3(XMIN_BOX, YMIN_BOX, ZMIN_BOX), // Point B
		glm::vec3(XMAX_BOX, YMIN_BOX, ZMIN_BOX), // Point C
		glm::vec3(XMAX_BOX, YMAX_BOX, ZMIN_BOX), // Point D
		glm::vec3(XMIN_BOX, YMAX_BOX, ZMIN_BOX)	 // Point A
	);
	frontWallPlane->setColor(glm::vec3(1, 1, 1)); // Set colour
	frontWallPlane->setReflectivity(false);
	sceneObjects.push_back(frontWallPlane);

	float CDR = 3.14159265 / 180.0; // Conversion from degrees to radians
	glm::mat4 transform = glm::mat4(1.0);
	transform = glm::translate(transform, glm::vec3(0.0, 0.05, 0.0));
	transform = glm::rotate(transform, 15 * CDR, glm::vec3(0.0, 0.0, -1.0));
	transform = glm::scale(transform, glm::vec3(0.5, 1, 0.5));

	Sphere *smallSphere = new Sphere(glm::vec3(0.0, 0.0, 0.0), 10.0);
	smallSphere->setColor(glm::vec3(1, 1, 1)); // Set colour to black
	// smallSphere->setReflectivity(true, 0.8);
	smallSphere->setTransform(true, transform);

	sceneObjects.push_back(smallSphere);

	Sphere *smallSphere2 = new Sphere(glm::vec3(20.0, -40, -30.0), 10.0);
	smallSphere2->setColor(glm::vec3(1, 1, 1)); // Set colour to black
	smallSphere2->setTransparency(true, 0.7);
	smallSphere2->setReflectivity(true, 0.1);
	sceneObjects.push_back(smallSphere2);

	Sphere *smallSphere3 = new Sphere(glm::vec3(-20.0, -40, -30.0), 10.0);
	smallSphere3->setColor(glm::vec3(1, 1, 1)); // Set colour to black
	smallSphere3->setRefractivity(true, 0.8, 0.9);
	sceneObjects.push_back(smallSphere3);

	Cylinder *cyl = new Cylinder(glm::vec3(30.0, -50, 30.0), 10.0, 20);
	cyl->setColor(glm::vec3(0, 0, 1)); // Set colour to blue
	sceneObjects.push_back(cyl);
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(1080, 1080);
	glutInitWindowPosition(20, 20);
	glutCreateWindow("Raytracing");

	glutDisplayFunc(display);
	initialize();

	glutMainLoop();
	return 0;
}
