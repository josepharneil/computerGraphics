#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits>
#include <math.h>
#include <random>
#include <time.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

#define SCREEN_WIDTH 400
#define SCREEN_HEIGHT 400
#define FULLSCREEN_MODE false
#define PI 3.14159265
#define RAYDEPTH 4
#define DIFFUSE_SAMPLES 1000000000
#define LIGHT_POWER 8.0f
// #define FOCAL_SPHERE_RADIUS 250.0f
#define APERTURE 0.2f
// #define isAAOn false

//============= Global Variables =============//
int t;
bool quit;


//============= Overrides =============//
#pragma region Overrides
//Print vec4s
std::ostream &operator<<( std::ostream &os, vec4 const &v )
{
  return os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w
	          << ")";
}

std::ostream &operator<<( std::ostream &os, vec3 const &v )
{
  return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

std::ostream &operator<<( std::ostream &os, mat4 const &m )
{
	glm::mat4 mt = transpose( m );
  return os << mt[ 0 ] << endl
	          << mt[ 1 ] << endl
	          << mt[ 2 ] << endl
	          << mt[ 3 ];
}
#pragma endregion


//============= Structures =============//
#pragma region Structures
struct Intersection
{
    vec4 position;
    float distance;
    int triangleIndex;
};

// struct Triangle
// {
//   vec4 v0;
//   vec4 v1;
//   vec4 v2;
//   vec4 normal;
//   vec3 color;
// };
#pragma endregion Structures

//============= Function Defs =============//
#pragma region FunctionDefs
void Update(vec4& cameraPos, int& yaw, vec4& lightPos, mat4& cameraMatrix, bool& isAAOn, 
                              vec3 screenAcc[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, 
                           int& yaw, vec4& lightPos, vec3& lightColour, mat4& cameraMatrix,bool& isAAOn, 
                           vec3 screenAccumulator[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad);
bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection);
vec3 DirectLight( Intersection& intersection, vec4& lightPos, 
                  vec3& lightColour, vector<Triangle>& triangles);
vec3 PathTracer(Intersection current, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles, int depth, vec3 previous, bool isSampleDirectLight );
void CreateCoordinateSystem(const vec3& N, vec3& Nt, vec3& Nb);
vec3 UniformSampleHemisphere(const float &rand1, const float &rand2);
vec3 Reflect(const vec3 &incident, const vec3 &normal);
void ResetScreenAccumulator(vec3 screenAcc[SCREEN_WIDTH][SCREEN_HEIGHT]);
vec3 Vec4ToVec3(vec4& vec4);
vec4 Vec3ToHomogenous(vec3& vec3);
void PrintPairOfNumbers(float f1, float f2);

#pragma endregion FunctionDefs


//============= Main =============//
int main( int argc, char* argv[] )
{
  //Initially, do not quit
  quit = false;

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  //Instantiate vector of triangles
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  vector<Triangle> originalTriangles;
  LoadTestModel( originalTriangles );

  //Camera control
  vec4 cameraPos(0.0f,0.0f,-1.8f,1.0f);
  // vec4 cameraPos(0.0f,0.0f,-1.8f,1.0f);
  mat4 cameraMatrix;
  int yaw = 0;

  //Focusing
  float focalSphereRad = 1.4f;

  //Initialise camera matrix
  cameraMatrix[0][0] = cos( yaw * PI / 180 );
  cameraMatrix[0][1] = 0;
  cameraMatrix[0][2] = sin( yaw * PI / 180 );
  cameraMatrix[0][3] = 0;

  cameraMatrix[1][0] = 0;
  cameraMatrix[1][1] = 1;
  cameraMatrix[1][2] = 0;
  cameraMatrix[1][3] = 0;

  cameraMatrix[2][0] = -sin( yaw * PI / 180 );
  cameraMatrix[2][1] = 0;
  cameraMatrix[2][2] = cos( yaw * PI / 180 );
  cameraMatrix[2][3] = 0;

  cameraMatrix[3][0] = cameraPos.x;
  cameraMatrix[3][1] = cameraPos.y;
  cameraMatrix[3][2] = cameraPos.z;
  cameraMatrix[3][3] = 1;
  // int yaw = 20;

  //Create light source
  vec4 lightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec4 originalLightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec3 lightColour = LIGHT_POWER * vec3( 1.0f, 1.0f, 1.0f );

  bool isAAOn = false;

  //Set up array of pixel vals
  vec3 screenAccumulator[SCREEN_WIDTH][SCREEN_HEIGHT];
  int sampleCount = 1;

  //Set all to 0
  ResetScreenAccumulator(screenAccumulator);

  //Update and draw
  while( !quit ) //NoQuitMessageSDL() )
  {
    Update(cameraPos, yaw, originalLightPos, cameraMatrix, isAAOn, screenAccumulator, sampleCount, focalSphereRad);

    //Rotation
    mat4 invCameraMatrix = glm::inverse(cameraMatrix);
    
    for (int t = 0; t < triangles.size(); t++)
    {
      triangles[t].v0 = invCameraMatrix * originalTriangles[t].v0;
      triangles[t].v1 = invCameraMatrix * originalTriangles[t].v1;
      triangles[t].v2 = invCameraMatrix * originalTriangles[t].v2;
      triangles[t].ComputeNormal();
    }

    lightPos = invCameraMatrix * originalLightPos;

    
    Draw(screen, triangles, cameraPos, yaw, lightPos, lightColour, cameraMatrix, isAAOn, screenAccumulator, sampleCount, focalSphereRad);


    SDL_Renderframe(screen);
  }

  //Output
  SDL_SaveImage( screen, "screenshot.bmp" );

  //Kill screen
  KillSDL(screen);
  return 0;

}

//============= Draw =============//
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, 
                          int& yaw, vec4& lightPos, vec3& lightColour, mat4& cameraMatrix, bool& isAAOn, 
                          vec3 screenAccumulator[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad)
{
  if(sampleCount == DIFFUSE_SAMPLES) {return;}
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // FOCAL LENGTH
  float focalLength = SCREEN_WIDTH/2;

  //Indirect lighting approximation
  // vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );

  //Instantiate closest intersection
  Intersection closestIntersection;

  // #pragma omp parallel for collapse(2)
  //For each pixel
  for (int row = 0; row < SCREEN_HEIGHT; row++)
  {
    for (int col = 0; col < SCREEN_WIDTH; col++)
    {
      if(isAAOn)
      {
        //============= Anti-Aliasing =============//

        //Array of 4 closests intersections
        Intersection closestIntersections[4];
        //Array of whether intersection has occurred
        bool isIntersections[4];
        //Array of subpixel directions
        vec4 directions[4];

        //Initialise random
        int colRand;
        int rowRand;

        float colRandF;
        float rowRandF;

        //Imagine screen is flattened to 1D array
        //Seed is the pixel number
        // int seed = (row * SCREEN_WIDTH) + col;

        // srand(seed);
        
        //Get random numbers
        colRand = rand() % 1001;
        rowRand = rand() % 1001;

        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //Top-left subpixel
        directions[0] = normalize(vec4(
          col-SCREEN_WIDTH/2 - (colRandF),
          row-SCREEN_HEIGHT/2 - (rowRandF),
          focalLength,
          1.0f));

        colRand = rand() % 1001;
        rowRand = rand() % 1001;

        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //Top-right subpixel
        directions[1] = normalize(vec4(
          col-SCREEN_WIDTH/2 - (colRandF),
          row-SCREEN_HEIGHT/2 + (rowRandF),
          focalLength,
          1.0f));

        colRand = rand() % 1001;
        rowRand = rand() % 1001;

        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //bottom-left subpixel
        directions[2] = normalize(vec4(
          col-SCREEN_WIDTH/2 + (colRandF),
          row-SCREEN_HEIGHT/2 - (rowRandF),
          focalLength,
          1.0f));

        colRand = rand() % 1001;
        rowRand = rand() % 1001;

        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //bottom-right subpixel
        directions[3] = normalize(vec4(
          col-SCREEN_WIDTH/2 + (colRandF),
          row-SCREEN_HEIGHT/2 + (rowRandF),
          focalLength,
          1.0f));


        //For each sub-pixel
        for (int i = 0; i < 4; i++)
        {
          //Compute ClosestIntersection
          isIntersections[i] = ClosestIntersection(
            vec4(0,0,0,1),
            directions[i],
            triangles,
            closestIntersections[i]);
        }

        vec3 colourTotal = vec3(0.0f,0.0f,0.0f);

        //For each subpixel
        for (int i = 0; i < 4; i++)
        {
          //If there is an intersection
          if (isIntersections[i])
          {
            //Get triangle from triangles
            Triangle intersectedTriangle = triangles[closestIntersections[i].triangleIndex];

            //Only sample direct light for 0th sample
            bool isSampleDirectLight;
            if(sampleCount == 1)
            {
              isSampleDirectLight = true;
            }
            else
            {
              isSampleDirectLight = false;
            }

            vec3 pathTracedLight = PathTracer(closestIntersections[i], lightPos, lightColour, triangles, 0, vec3(0,0,0), isSampleDirectLight);

            // vec3 colour = pathTracedLight;// * intersectedTriangle.color;

            //Running colour total over whole pixel
            colourTotal += pathTracedLight;
          }
        }
        //average colour over 4 sub-pixels
        colourTotal = colourTotal/4.0f;
        // colourTotal = vec3(colourTotal.x/4.0f,colourTotal.y/4.0f,colourTotal.z/4.0f)
        // cout << colourTotal << "\n";

        //Accumulate
        screenAccumulator[col][row] += colourTotal;

        //Average over #samples
        vec3 currentColour = vec3(screenAccumulator[col][row].x/sampleCount,screenAccumulator[col][row].y/sampleCount,screenAccumulator[col][row].z/sampleCount);

        //set to colour of that triangle
        PutPixelSDL(screen, col, row, currentColour);
      }

      else//No AA
      {
        //Compute ray direction
        vec4 direction = vec4(
          col-SCREEN_WIDTH/2,
          row-SCREEN_HEIGHT/2,
          focalLength,
          1.0f);

        //Normalise direction of ray
        direction = normalize(direction);


        vec4 focalPoint = direction * focalSphereRad;
        // cout << focalSphereRad << "\n";
        focalPoint.w = 1.0f;

        float randX = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);
        float randY = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);
        float randZ = 0.0f;//APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);

        direction = normalize(focalPoint - vec4(randX,randY,randZ,0.0f));

        //Compute ClosestIntersection
        bool intersect = ClosestIntersection(
          vec4(randX,randY,randZ,1.0f),
          direction,
          triangles,
          closestIntersection);

        //If an intersection occurs
        if (intersect)
        {
          //Get triangle from triangles
          Triangle intersectedTriangle = triangles[closestIntersection.triangleIndex];

          //Only sample direct light for 0th sample
          bool isSampleDirectLight;
          if(sampleCount == 1)
          {
            isSampleDirectLight = true;
          }
          else
          {
            isSampleDirectLight = false;
          }
          
          //Path trace the light
          vec3 pathTracedLight = PathTracer(closestIntersection, lightPos, lightColour, triangles, 0, vec3(0,0,0), isSampleDirectLight);

          //Accumulate
          screenAccumulator[col][row] += pathTracedLight;

          //Average over #samples
          vec3 currentColour = vec3(screenAccumulator[col][row].x/sampleCount,screenAccumulator[col][row].y/sampleCount,screenAccumulator[col][row].z/sampleCount);
          //set to colour of that triangle
          PutPixelSDL(screen, col, row, currentColour);
        }
      }//end no aa
    }
  }//end loop through pixels

  //Increment sample
  sampleCount += 1;
}

//============= Update =============//
/*Place updates of parameters here*/
void Update(vec4& cameraPos, int& yaw, vec4& lightPos, mat4& cameraMatrix, bool& isAAOn, 
                              vec3 screenAcc[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad)
{
  // static int t = SDL_GetTicks();
  /* Compute frame time */
  // int t2 = SDL_GetTicks();
  // float dt = float(t2-t);
  // t = t2;

  float movementSpeed = 0.1f;

  SDL_Event e;

  while(SDL_PollEvent(&e))
  {
    if( e.type != SDL_KEYDOWN) {continue;}

    //Move camera with arrow keys

    //MOVE FORWARDS/BACKWARDS
    if( e.key.keysym.scancode == SDL_SCANCODE_UP )
    {
      cameraPos.z += movementSpeed;

    }
    if( e.key.keysym.scancode == SDL_SCANCODE_DOWN )
    {
      //Move camera backward
      cameraPos.z -= movementSpeed;

    }

    //MOVE LEFT/RIGHT
    if( e.key.keysym.scancode == SDL_SCANCODE_LEFT )
    {
      //Move camera left
      cameraPos.x -= movementSpeed;
    }
    if( e.key.keysym.scancode == SDL_SCANCODE_RIGHT )
    {
      //Move camera right
      cameraPos.x += movementSpeed;
    }

    //MOVE UP/DOWN
    if( e.key.keysym.scancode == SDL_SCANCODE_J )
    {
      //Move camera left
      cameraPos.y -= movementSpeed;
    }
    if( e.key.keysym.scancode == SDL_SCANCODE_M )
    {
      //Move camera right
      cameraPos.y += movementSpeed;
    }

    //Move light with WASD
    if ( e.key.keysym.scancode == SDL_SCANCODE_W )
    {
      lightPos.z += movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_S )
    {
      lightPos.z -= movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_D )
    {
      lightPos.x += movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_A )
    {
      lightPos.x -= movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_R )
    {
      lightPos.y -= movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_F )
    {
      lightPos.y += movementSpeed;
    }

    //ROTATE AROUND Y-AXIS
    if ( e.key.keysym.scancode == SDL_SCANCODE_PERIOD )
    {
      yaw -= 5;
      yaw = yaw % 360;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_COMMA )
    {
      yaw += 5;
      yaw = yaw % 360;
    }

    //Switch on/off AA
    if ( e.key.keysym.scancode == SDL_SCANCODE_N )
    {
      if(isAAOn)
      {
        isAAOn = false;
      }
      else
      {
        isAAOn = true;
      }
    }

    if (e.key.keysym.scancode == SDL_SCANCODE_O)
    {
      focalSphereRad += 0.2f;
    }
      if (e.key.keysym.scancode == SDL_SCANCODE_L)
    {
      focalSphereRad -= 0.2f;
    }

    //Restart sampling on any input
    ResetScreenAccumulator(screenAcc);
    sampleCount = 1;

    //Recompute cameraMat on any input
    cameraMatrix[0][0] = cos( yaw * PI / 180 );
    cameraMatrix[0][1] = 0;
    cameraMatrix[0][2] = sin( yaw * PI / 180 );
    cameraMatrix[0][3] = 0;

    cameraMatrix[1][0] = 0;
    cameraMatrix[1][1] = 1;
    cameraMatrix[1][2] = 0;
    cameraMatrix[1][3] = 0;

    cameraMatrix[2][0] = -sin( yaw * PI / 180 );
    cameraMatrix[2][1] = 0;
    cameraMatrix[2][2] = cos( yaw * PI / 180 );
    cameraMatrix[2][3] = 0;

    cameraMatrix[3][0] = cameraPos.x;
    cameraMatrix[3][1] = cameraPos.y;
    cameraMatrix[3][2] = cameraPos.z;
    cameraMatrix[3][3] = 1;


    //Quit trigger
    if( e.type == SDL_QUIT )
    {
      quit = true;
    }
    if( e.type == SDL_KEYDOWN )
    {
      if( e.key.keysym.sym == SDLK_ESCAPE)
      {
        quit = true;
      }
    }
  }
}


//======== Raycasting ========//
bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection)
{
  start = start + 1e-4f * dir;
  //Get closest intersection
  float minDist = numeric_limits<float>::max();
  bool result = false;

  //For each triangle
  for (int i = 0; i < triangles.size(); i++)
  {
    Triangle triangle = triangles[i];

    //3D version of dir (ignore last homogenous component)
    vec3 dir3 = vec3(dir.x, dir.y, dir.z);

    //Define vertices
    vec4 v0 = triangle.v0;
    vec4 v1 = triangle.v1;
    vec4 v2 = triangle.v2;

    //Define axes
    vec3 e1 = vec3(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z);

    //b = s - v0
    vec3 b = vec3(start.x - v0.x, start.y - v0.y, start.z - v0.z);

    //A
    mat3 A = mat3( -dir3, e1, e2 );

    //x = (t u v)^T
    //COME BACK HERE todo cramar's rule instead of inbuilt inverse
    vec3 x = glm::inverse( A ) * b;

    float t = x.x;//if neg, can reject immediately
    float u = x.y;
    float v = x.z;

    //Check if is intersection is within triangle
    if ( (u >= 0.0f) && (v >= 0.0f) && (u + v <= 1.001f) && (t > 0.0f) )
    {
      result = true;
      //Compute "intersection" structure
      if (t < minDist)
      {
        //We set w = 1 -- not sure of this
        closestIntersection.position = vec4(start.x+t*dir.x, start.y+t*dir.y, start.z+t*dir.z, 1);
        closestIntersection.distance = t;
        closestIntersection.triangleIndex = i;

        minDist = t;
      }
    }
  }
  return result;
}



//======== Lighting ========//
//Returns power only - colour is ignored
vec3 DirectLight( Intersection& intersection, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles )
{
  // Light colour is P
  vec3 P = lightColour;
  // Get normal to triangle
  vec4 nNorm = normalize(triangles[intersection.triangleIndex].normal);
  // r is vector from intersection point to light source
  vec4 r  = lightPos - intersection.position;
  vec4 rNorm = normalize(r);

  //Compute power per area B
  float A = 4.0f * M_PI * ( pow(glm::length(r),2.0f) );

  vec3 B = vec3(P.x/A,P.y/A,P.z/A);
 
  //Compute power per real surface D
  vec3 D = B * max(glm::dot (rNorm,nNorm) , 0.0f);

  //Vector direction from surface to light
  vec4 lightDir = lightPos - intersection.position;

  vec4 lightDirNorm = normalize(lightPos - intersection.position);

  //Distance from surface to light
  float lightDist = glm::length(lightDir);

  //Intersection of ray cast from surface point to light
  Intersection lightIntersection;

  //Compute light intersection
  //Move all points a very small amount along the ray.
  //This 1e-5f is to account for these floating point errors.
  // bool intersectionFound = ClosestIntersection(intersection.position + 1e-4f * lightDirNorm,lightDirNorm,triangles,lightIntersection);
  bool intersectionFound = ClosestIntersection(intersection.position,lightDirNorm,triangles,lightIntersection);
  
  //If an intersection is found
  //2nd part of conjunction is to deal with the light being inside of the wall 
  //(when the light rays are parallel to the wall, i.e., perpendicular to the normal))
  if (intersectionFound && ( abs(glm::dot( nNorm, (intersection.position - lightPos))) >= 0 ) )
  {
    //If the ray is stopped short by another intersection
    if ( glm::length(lightIntersection.position - intersection.position) < lightDist )
    {
      //Set to black (i.e., a shadow)
      D = vec3(0.0f,0.0f,0.0f);
    }
  }

  return D;//(D * triangles[intersection.triangleIndex].color);
}

//current is the point we want to find lighting for; previous is where we cast a ray from to find this point.
vec3 PathTracer(Intersection current, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles, int depth, vec3 previous, bool isSampleDirectLight )
{
  // int chamberSize = 10;
  // int russianRoulette = rand() % (chamberSize + 1);
  // if(russianRoulette == chamberSize && depth > 1)
  // {
  //   return vec3(0,0,0);
  // }

  //Stop recursion if depth exceed
  if(depth > RAYDEPTH) 
  { 
    return vec3(0,0,0);
  }

  depth += 1;
 
  //Initialise result to be returned
  vec3 result = vec3(0.0f,0.0f,0.0f);

  //============= Smoothness =============//
  //If there is any smoothness (for reflectance)
  if(triangles[current.triangleIndex].smoothness == 1.0f)
  {
    //Find reflected ray from indidence ray and normal
    vec3 incidentRay = Vec4ToVec3(current.position) - previous;
    vec3 normalVec3 = Vec4ToVec3(triangles[current.triangleIndex].normal);
    vec3 reflectedRay = Reflect(incidentRay, normalVec3);

    //Find the next intersection of the reflected ray
    Intersection nextIntersection;
    vec4 reflectedRayV4 = Vec3ToHomogenous(reflectedRay);
    bool isIntersection = ClosestIntersection(current.position,reflectedRayV4,triangles,nextIntersection);

    //If there is an intersection
    if(isIntersection)
    {
      //Add (slight attenuated) directlight from subsequent reflections, taking mirror colour into account
      result += (0.8f * PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,incidentRay,true) * triangles[current.triangleIndex].color);
    }
    else//If there is no intersection, don't recurse anymore 
    {
      //(effectively the reflected ray shoots into empty space)
      // result += vec3(0.0f,0.0f,0.0f);
    }

    //If we want the mirror itself to have shadows, colour etc.
    //Then backtrace light from the mirror surface
    // result += (DirectLight( current, lightPos, lightColour, triangles ) * triangles[current.triangleIndex].color);
    
  }
  //============= End Smoothness =============//


  //============= Diffuse =============//
  //If the surface is diffuse (no smoothness)
  if(triangles[current.triangleIndex].smoothness == 0.0f)
  {
    //indirectlight comes from sampling
    vec3 indirectLight;

    //Find local coordinate system at the intersection
    vec3 normal = Vec4ToVec3(triangles[current.triangleIndex].normal);
    vec3 Nt;
    vec3 Nb;
    CreateCoordinateSystem(normal,Nt,Nb);

    float PDF = 1.0f/(2.0f*PI);


    //Sample two random numbers from 0->1
    // srand(time( NULL ));
    float rand1 = ((float) rand() / (RAND_MAX));
    float rand2 = ((float) rand() / (RAND_MAX));

    //Local hemispherical sampling
    vec3 LocalSampleDir = UniformSampleHemisphere(rand1,rand2);

    //Transform to global
    vec4 WorldSampleDir = vec4(LocalSampleDir.x * Nb.x + LocalSampleDir.y * normal.x + LocalSampleDir.z * Nt.x, 
                               LocalSampleDir.x * Nb.y + LocalSampleDir.y * normal.y + LocalSampleDir.z * Nt.y, 
                               LocalSampleDir.x * Nb.z + LocalSampleDir.y * normal.z + LocalSampleDir.z * Nt.z,1.0f);
    
    //"Cast" ray, ie., find ClosestIntersection
    Intersection nextIntersection;
    bool isIntersect = ClosestIntersection(current.position,WorldSampleDir,triangles,nextIntersection);

    //If there is an interesction, recurse at +1 depth
    if(isIntersect)
    {
      // float cosTheta = 2.0f * PI * rand2;
      // indirectLight += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position)) * (cosTheta));

      float cosTheta = rand1;
      indirectLight += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true)*cosTheta*0.5f);
    }
    else//If there is no intersection, do nothing
    {
      //do nothing
    }

    //Div. indirect light by pdf
    indirectLight = (indirectLight/PDF);


    //Backtrace light from the (non-reflective) surface
    result += (indirectLight) * triangles[current.triangleIndex].color;
  }
  //============= End Diffuse =============//



  //============= Emission =============//
  result += triangles[current.triangleIndex].emissive;
  //============= End Emission =============//


  //Direct lighting
  // if(isSampleDirectLight)
  if(true)
  {
    // float temp = 2.0f * PI;
    result += ( (DirectLight( current, lightPos, lightColour, triangles ) * triangles[current.triangleIndex].color));
  }


  //Return result
  return result;
}


vec3 Reflect(const vec3 &incident, const vec3 &normal)
{
  return (incident - ((2.0f * glm::dot(incident,normal)) * normal));
}


void CreateCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb) 
{ 
    if (fabs(N.x) > fabs(N.y)) 
    {
      Nt = vec3(N.z, 0.0f, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
    }
    else
    {
      Nt = vec3(0.0f, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z); 
    }
    Nb = glm::cross(N,Nt);
}

vec3 UniformSampleHemisphere(const float &rand1, const float &rand2) 
{ 
    // cos(theta) = r1 = y
    // cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
    float sinTheta = sqrtf(1.0f - rand1 * rand1);
    float phi = 2 * M_PI * rand2; 
    float x = sinTheta * cosf(phi); 
    float z = sinTheta * sinf(phi); 
    return vec3(x, rand1, z);


  // const float r = sqrtf(rand1);
  // const float theta = 2.0f * PI * rand2;
  //   const float x = r * cosf(theta);
  //   const float y = r * sinf(theta);
 
  //   return vec3(x, y, sqrtf(max(0.0f, 1.0f - rand1)));
} 


void ResetScreenAccumulator(vec3 screenAcc[SCREEN_WIDTH][SCREEN_HEIGHT])
{
  //Set all to 0
  for(int i = 0; i < SCREEN_WIDTH; i++)
  {
    for(int j = 0; j < SCREEN_HEIGHT; j++)
    {
      screenAcc[i][j] = vec3(0,0,0);
    }
  }
}



//Helpful functions
vec3 Vec4ToVec3(vec4& vec4)
{
  return vec3(vec4.x,vec4.y,vec4.z);
}

vec4 Vec3ToHomogenous(vec3& vec3)
{
  return vec4(vec3.x,vec3.y,vec3.z,1);
}

void PrintPairOfNumbers(float f1, float f2)
{
  cout << "(" << f1 << "," << f2 << ")\n";
}
