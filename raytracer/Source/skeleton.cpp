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
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

#define SCREEN_WIDTH 400
#define SCREEN_HEIGHT 400
#define FULLSCREEN_MODE false
#define PI 3.14159265
#define RAYDEPTH 10
#define DIFFUSE_SAMPLES 1000000000
#define LIGHT_POWER 4.5f
#define SPHERE_LIGHT_RADIUS 0.05f
#define SPHERE_LIGHT_SAMPLES 10
// #define FOCAL_SPHERE_RADIUS 250.0f
#define APERTURE 0.0f//e.g. 0.1f
// #define isAAOn false
#define FOG_STRENGTH 0.1f

//============= Global Variables =============//
int timeT;
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
  int sphereIndex;
};

// struct Triangle
// {
//   vec4 v0;
//   vec4 v1;
//   vec4 v2;
//   vec4 normal;
//   vec3 color;
// };

// struct Sphere
// {
//   vec4 centre;
//   float radius;
//  vec3 color;
//  vec3 emissive;
//  float smoothness;
// };

class Sphere
{
  public:
    vec4 centre;
    float radius;
    vec3 color;
    vec3 emissive;
    float smoothness;
    float refractiveIndex;

  Sphere( vec4 centre, float radius, vec3 color, vec3 emissive, float smoothness, float refractiveIndex )
    : centre(centre), radius(radius), color(color), emissive(emissive), smoothness(smoothness), refractiveIndex(refractiveIndex){}
};




#pragma endregion Structures

//============= Function Defs =============//
#pragma region FunctionDefs
void Update(vec4& cameraPos, int& yaw, vec4& lightPos, mat4& cameraMatrix, bool& isAAOn, 
                              vec3 screenAcc[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad, bool& isAreaLight,bool& isFog);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, 
                           int& yaw, vec4& lightPos, vec3& lightColour, mat4& cameraMatrix,bool& isAAOn, 
                           vec3 screenAccumulator[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad, bool& isAreaLight,
                           const vector<Sphere>& spheres,bool& isFog);
bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection, const vector<Sphere>& spheres);
vec3 DirectLight( Intersection& intersection, vec4& lightPos, 
                  vec3& lightColour, vector<Triangle>& triangles, const vector<Sphere>& spheres);
vec3 AreaLightSample( Intersection& intersection, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles,const vector<Sphere>& spheres );
vec3 PathTracer(Intersection current, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles, int depth, vec3 previous, bool isSampleDirectLight, bool &isAreaLight,const vector<Sphere>& spheres );
void CreateCoordinateSystem(const vec3& N, vec3& Nt, vec3& Nb);
vec3 UniformSampleHemisphere(const float &rand1, const float &rand2);
float FogAmount(const vec3 start,const vec3 end);
float FogUnitCurveArea(float end, float start, float amplitude, float frequency);
vec3 Refract(const vec3 &incidentRay, const vec3 &normal, const float &IoR);
vec3 Reflect(const vec3 &incident, const vec3 &normal);
void ResetScreenAccumulator(vec3 screenAcc[SCREEN_WIDTH][SCREEN_HEIGHT]);
vec4 NormaliseNoHomogenous(vec4 vector4);
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
  timeT = SDL_GetTicks(); /*Set start value for timer.*/

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

  bool isAreaLight = false;
  
  bool isFog = false;

  

  //Set up array of pixel vals
  vec3 screenAccumulator[SCREEN_WIDTH][SCREEN_HEIGHT];
  int sampleCount = 1;

  //Set all to 0
  ResetScreenAccumulator(screenAccumulator);


  //============= SPHERES =============//
  vector<Sphere> spheres;
  vector<Sphere> originalSpheres;

  spheres.clear();
  spheres.reserve( 5*2*3 );

  originalSpheres.clear();
  originalSpheres.reserve( 5*2*3 );

  vec4 centreTEMP(-0.4f,-0.1f,-0.5f,1.0f);
  vec4 centreTEMP2(0.3f,-0.2f,-0.3f,1.0f);
  Sphere mySphere = Sphere(centreTEMP,0.2f,vec3(0.75f,0.75f,0.75f), vec3(0.0f,0.0f,0.0f),0.7f,1.0f);
  Sphere mySphere2 = Sphere(centreTEMP2,0.2f,vec3(1.0f,1.0f,1.0f), vec3(0.0f,0.0f,0.0f),0.0f,2.42f);
  spheres.push_back( mySphere );
  spheres.push_back( mySphere2 );
  originalSpheres.push_back( mySphere );
  originalSpheres.push_back( mySphere2 );


  //Update and draw
  while( !quit ) //NoQuitMessageSDL() )
  {
    Update(cameraPos, yaw, originalLightPos, cameraMatrix, isAAOn, screenAccumulator, sampleCount, focalSphereRad, isAreaLight,isFog);

    //Rotation
    mat4 invCameraMatrix = glm::inverse(cameraMatrix);
    
    for (int t = 0; t < triangles.size(); t++)
    {
      // cout << triangles[t].v0 <<"\n";
      triangles[t].v0 = invCameraMatrix * originalTriangles[t].v0;
      triangles[t].v1 = invCameraMatrix * originalTriangles[t].v1;
      triangles[t].v2 = invCameraMatrix * originalTriangles[t].v2;
      triangles[t].ComputeNormal();
    }

    for (int s = 0; s < spheres.size(); s++)
    {
      spheres[s].centre = invCameraMatrix * originalSpheres[s].centre;
    }

    // cout << spheres[0].centre << "\n";

    lightPos = invCameraMatrix * originalLightPos;

    
    Draw(screen, triangles, cameraPos, yaw, lightPos, lightColour, cameraMatrix, isAAOn, screenAccumulator, sampleCount, focalSphereRad, isAreaLight, spheres,isFog);


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
                          vec3 screenAccumulator[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad, bool& isAreaLight, const vector<Sphere>& spheres, bool& isFog)
{
  if(sampleCount == DIFFUSE_SAMPLES) {return;}
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // FOCAL LENGTH
  float focalLength = SCREEN_WIDTH/1.5f;
  cout << focalLength << "\n";

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
        //Random aperture sample
        vec2 apertureSample[4];
        float randX;
        float randY;
        vec4 focalPoint;

        //Initialise random
        int colRand;
        int rowRand;

        float colRandF;
        float rowRandF;

        //Imagine screen is flattened to 1D array
        //Seed is the pixel number
        // int seed = (row * SCREEN_WIDTH) + col;

        // srand(seed);
        

        /////Subpixel1
        colRand = rand() % 1001;
        rowRand = rand() % 1001;
        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //Top-left subpixel
        directions[0] = NormaliseNoHomogenous(vec4(
          col-SCREEN_WIDTH/2 - (colRandF),
          row-SCREEN_HEIGHT/2 - (rowRandF),
          focalLength,
          1.0f));

        directions[0] = NormaliseNoHomogenous(directions[0]);

        focalPoint = directions[0] * focalSphereRad;
        focalPoint.w = 1.0f;

        randX = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);
        randY = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);

        directions[0] = NormaliseNoHomogenous(focalPoint - vec4(randX,randY,0.0,0.0f));

        apertureSample[0].x = randX;
        apertureSample[0].y = randY;




        /////Subpixel2
        colRand = rand() % 1001;
        rowRand = rand() % 1001;
        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //Top-right subpixel
        directions[1] = NormaliseNoHomogenous(vec4(
          col-SCREEN_WIDTH/2 - (colRandF),
          row-SCREEN_HEIGHT/2 + (rowRandF),
          focalLength,
          1.0f));

        directions[1] = NormaliseNoHomogenous(directions[1]);

        focalPoint = directions[1] * focalSphereRad;
        focalPoint.w = 1.0f;

        randX = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);
        randY = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);

        directions[1] = NormaliseNoHomogenous(focalPoint - vec4(randX,randY,0.0,0.0f));

        apertureSample[1].x = randX;
        apertureSample[1].y = randY;


        /////Subpixel3
        colRand = rand() % 1001;
        rowRand = rand() % 1001;
        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //bottom-left subpixel
        directions[2] = NormaliseNoHomogenous(vec4(
          col-SCREEN_WIDTH/2 + (colRandF),
          row-SCREEN_HEIGHT/2 - (rowRandF),
          focalLength,
          1.0f));

        directions[2] = NormaliseNoHomogenous(directions[2]);

        focalPoint = directions[2] * focalSphereRad;
        focalPoint.w = 1.0f;

        randX = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);
        randY = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);

        directions[2] = NormaliseNoHomogenous(focalPoint - vec4(randX,randY,0.0,0.0f));

        apertureSample[2].x = randX;
        apertureSample[2].y = randY;


        /////Subpixel4
        colRand = rand() % 1001;
        rowRand = rand() % 1001;
        colRandF = float(colRand) / 2000.0f;
        rowRandF = float(rowRand) / 2000.0f;

        //bottom-right subpixel
        directions[3] = NormaliseNoHomogenous(vec4(
          col-SCREEN_WIDTH/2 + (colRandF),
          row-SCREEN_HEIGHT/2 + (rowRandF),
          focalLength,
          1.0f));

        directions[3] = NormaliseNoHomogenous(directions[3]);

        focalPoint = directions[3] * focalSphereRad;
        focalPoint.w = 1.0f;

        randX = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);
        randY = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);

        directions[3] = NormaliseNoHomogenous(focalPoint - vec4(randX,randY,0.0,0.0f));

        apertureSample[3].x = randX;
        apertureSample[3].y = randY;



        //For each sub-pixel
        for (int i = 0; i < 4; i++)
        {
          //Compute ClosestIntersection
          isIntersections[i] = ClosestIntersection(
            vec4(apertureSample[i].x,apertureSample[i].y,0.0f,1),
            directions[i],
            triangles,
            closestIntersections[i],spheres);
        }

        vec3 colourTotal = vec3(0.0f,0.0f,0.0f);

        //For each subpixel
        for (int i = 0; i < 4; i++)
        {
          //If there is an intersection
          if (isIntersections[i])
          {
            //Get triangle from triangles
            // Triangle intersectedTriangle = triangles[closestIntersections[i].triangleIndex];

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

            vec3 pathTracedLight = PathTracer(closestIntersections[i], lightPos, lightColour, triangles, 0, vec3(0,0,0), isSampleDirectLight,isAreaLight,spheres);

            if(isFog)
            {            
              //============= Fog =============//
              float fogAmount = FogAmount(vec3(apertureSample[i].x,apertureSample[i].y,0.0f),Vec4ToVec3(closestIntersections[i].position));
              vec3 fogColour = vec3(0.8f*fogAmount,0.8f*fogAmount,0.8f*fogAmount);
              // cout << fogAmount << "\n";
              colourTotal += fogColour;
              //============= END Fog =============//
            }

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
        direction = NormaliseNoHomogenous(direction);
        // if (/col == (SCREEN_WIDTH/2) && row == (SCREEN_HEIGHT/2)) {
          // cout << direction << endl;
        // }
        //Get focal length vector
        vec4 focalPoint = direction * focalSphereRad;
        focalPoint.w = 1.0f;

        //Randomly sample
        float randX = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);
        float randY = APERTURE * (((float) rand() / (RAND_MAX)) - 0.5f);

        //New direction
        direction = NormaliseNoHomogenous(focalPoint - vec4(randX,randY,0.0f,0.0f));


        //Compute ClosestIntersection
        bool intersect = ClosestIntersection(
          vec4(randX,randY,0.0f,1.0f),
          direction,
          triangles,
          closestIntersection,
          spheres);

        // cout << spheres.size() << "\n";
        

        //If an intersection occurs
        if (intersect)
        {
          // cout << spheres.size() << "\n";
          //Get triangle from triangles
          // Triangle intersectedTriangle = triangles[closestIntersection.triangleIndex];
          // if(closestIntersection.triangleIndex != -1)
          // {
          //   PutPixelSDL(screen, col, row, vec3(0.25f,0.0f,0.25f));
          // }
          // if(closestIntersection.sphereIndex != -1)
          // {
          //   PutPixelSDL(screen, col, row, vec3(1.0f,1.0f,1.0f));
          // }

          // Only sample direct light for 0th sample
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
          vec3 pathTracedLight = PathTracer(closestIntersection, lightPos, lightColour, triangles, 0, vec3(0,0,0), isSampleDirectLight,isAreaLight,spheres);

          //Accumulate
          screenAccumulator[col][row] += pathTracedLight;

          // cout << cameraPos << "\n";

          //Average over #samples
          vec3 currentColour = vec3(screenAccumulator[col][row].x/sampleCount,screenAccumulator[col][row].y/sampleCount,screenAccumulator[col][row].z/sampleCount);

          if(isFog)
          {
            //============= Fog =============//
            float fogAmount = FogAmount(vec3(randX,randY,0.0f),Vec4ToVec3(closestIntersection.position));
            vec3 fogColour = vec3(0.8f*fogAmount,0.8f*fogAmount,0.8f*fogAmount);
            // cout << fogAmount << "\n";
            currentColour += fogColour;
          //============= END Fog =============//
          }
          
          
          //set to colour of that triangle
          PutPixelSDL(screen, col, row, currentColour);
        }
      }//end no aa
    }
  }//end loop through pixels

  //Increment sample
  sampleCount += 1;
}

float FogAmount(const vec3 start,const vec3 end)
{
  float distance = glm::length(start - end);
  float baseAmount = distance * 0.1f;
  float varyingAmount = 0.0f;//initialise
  varyingAmount += FogUnitCurveArea(end.x,start.x,FOG_STRENGTH*0.6f,10.0f); //f was 0.6
  varyingAmount += FogUnitCurveArea(end.y,start.y,FOG_STRENGTH*1.2f,11.0f);   //f was 1.2 
  varyingAmount += FogUnitCurveArea(end.z,start.z,FOG_STRENGTH*0.9f,9.0f);   //f was 0.9
  varyingAmount *= distance;
  float resultAmount = baseAmount + varyingAmount;
  return min(resultAmount,1.0f);
}

float FogUnitCurveArea(float end,float start,float amplitude,float frequency)
{
  return (((end - (cos(frequency * end)/frequency)) * amplitude) - ((start - (cos(frequency * start)/ frequency)) * amplitude));
}

//============= Update =============//
/*Place updates of parameters here*/
void Update(vec4& cameraPos, int& yaw, vec4& lightPos, mat4& cameraMatrix, bool& isAAOn, 
                              vec3 screenAcc[SCREEN_WIDTH][SCREEN_HEIGHT], int &sampleCount, float& focalSphereRad, bool& isAreaLight, bool &isFog)
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
      // cout << cameraPos << "\n";

    }
    if( e.key.keysym.scancode == SDL_SCANCODE_DOWN )
    {
      //Move camera backward
      cameraPos.z -= movementSpeed;
      // cout << cameraPos << "\n";

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

    //Toggle arealight/pointlight
    if( e.key.keysym.scancode == SDL_SCANCODE_G)
    {
      isAreaLight = !isAreaLight;
    }

    //Move focalsphereradius
    if (e.key.keysym.scancode == SDL_SCANCODE_O)
    {
      focalSphereRad += 0.2f;
    }
    if (e.key.keysym.scancode == SDL_SCANCODE_L)
    {
      focalSphereRad -= 0.2f;
    }

    if (e.key.keysym.scancode == SDL_SCANCODE_L)
    {
      focalSphereRad -= 0.2f;
    }

    if (e.key.keysym.scancode == SDL_SCANCODE_Q)
    {
      isFog = !isFog;
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
  Intersection& closestIntersection,
  const vector<Sphere>& spheres)
{
  start = start + 1e-4f * dir;
  //Get closest intersection
  float minDist = numeric_limits<float>::max();
  bool result = false;

  //============= Triangles =============//
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
        closestIntersection.sphereIndex = -1;
        closestIntersection.triangleIndex = i;

        minDist = t;
      }
    }
    // return result;
  }
  //cout << dir << endl;
  //============= Spheres =============//
  //For each sphere
  for(int i = 0; i < spheres.size(); i++)
  {
    // cout << spheres[i].centre;
    float t0,t1;
    vec4 L = spheres[i].centre - start;
    // if (abs(dir.x) < 0.01f && abs(dir.y) < 0.01f)
      // cout << spheres[i].centre << " " << start << " " << dir << endl;
    float tca = glm::dot(L,dir);
    // cout << glm::dot(dir, dir) << "\n";
    if(tca < 0.0f)
    {
      continue;
    }
    float d2 = glm::dot(L,L) - (tca * tca); 
    
    // cout << d2<< "cont\n";
    if (d2 > (spheres[i].radius * spheres[i].radius))
    {
      continue;
    }
    //cout << (tca * tca) << "\n";
    float thc = sqrt((spheres[i].radius* spheres[i].radius) - d2); 
    t0 = tca - thc; 
    t1 = tca + thc; 

    if (t0 > t1) 
    {
      std::swap(t0, t1);
    }

    if (t0 < 0.0f)
    { 
      t0 = t1; // if t0 is negative, let's use t1 instead 
      if (t0 < 0.0f) 
      {
        continue;
      } // both t0 and t1 are negative 
    }

    //If get to this point then an intersection has occured.
    float t = t0;
    if (t < minDist)
    {
        //cout<<"sphere found" << "\n";
        closestIntersection.position = vec4(start.x+t*dir.x, start.y+t*dir.y, start.z+t*dir.z, 1);
        closestIntersection.distance = t;
        closestIntersection.sphereIndex = i;
        closestIntersection.triangleIndex = -1;

        // cout << closestIntersection.distance << "\n";
        minDist = t;
    }
    result = true;
  }

  return result;
}



//======== Lighting ========//
//Returns power only - colour is ignored
vec3 DirectLight( Intersection& intersection, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles,const vector<Sphere>& spheres )
{
  // Light colour is P
  vec3 P = lightColour;
  // Get normal to triangle
  vec4 nNorm;
  if(intersection.triangleIndex !=-1)
  {
    nNorm = NormaliseNoHomogenous(triangles[intersection.triangleIndex].normal);
    // return vec3(0,0,0);
  }
  else if(intersection.sphereIndex!= -1)
  {
    // cout <<"sphere\n";
    nNorm = NormaliseNoHomogenous(intersection.position - spheres[intersection.sphereIndex].centre);
  }

  
  // r is vector from intersection point to light source
  vec4 r  = lightPos - intersection.position;
  vec4 rNorm = NormaliseNoHomogenous(r);

  //Compute power per area B
  float A = 4.0f * M_PI * ( pow(glm::length(r),2.0f) );

  vec3 B = vec3(P.x/A,P.y/A,P.z/A);
 
  //Compute power per real surface D
  vec3 D = B * max(glm::dot (rNorm,nNorm) , 0.0f);

  //Vector direction from surface to light
  vec4 lightDir = lightPos - intersection.position;

  vec4 lightDirNorm = NormaliseNoHomogenous(lightPos - intersection.position);

  //Distance from surface to light
  float lightDist = glm::length(lightDir);

  //Intersection of ray cast from surface point to light
  Intersection lightIntersection;

  //Compute light intersection
  //Move all points a very small amount along the ray.
  //This 1e-5f is to account for these floating point errors.
  // bool intersectionFound = ClosestIntersection(intersection.position + 1e-4f * lightDirNorm,lightDirNorm,triangles,lightIntersection);
  bool intersectionFound = ClosestIntersection(intersection.position,lightDirNorm,triangles,lightIntersection,spheres);
  
  //If an intersection is found
  //2nd part of conjunction is to deal with the light being inside of the wall 
  //(when the light rays are parallel to the wall, i.e., perpendicular to the normal))
  if (intersectionFound && ( abs(glm::dot( nNorm, (intersection.position - lightPos))) >= 0 ) )
  {
    //If the ray is stopped short by another intersection
    if ( glm::length(lightIntersection.position - intersection.position) < lightDist )
    {
      //Set to black (i.e., a shadow)
      // D = vec3(1.0f,0.0f,1.0f);
    }
  }

  return D;//(D * triangles[intersection.triangleIndex].color);
}

//current is the point we want to find lighting for; previous is where we cast a ray from to find this point.
vec3 PathTracer(Intersection current, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles, int depth, vec3 previous, bool isSampleDirectLight, bool& isAreaLight, const vector<Sphere>& spheres )
{
  // float BRDF = 1.0f;
  //============= Russian Roulette =============//
  int chamberNumber = 20;
  int russianRoulette = rand() % (chamberNumber + 1);//number from 0->chamberNumber
  if(russianRoulette <= 6 && depth > 2)//60% chance, except at depth 0
  {
    return vec3(0,0,0);
  }
  //============= End Russian Roulette =============//

  //Stop recursion if depth exceed
  if(depth > RAYDEPTH) 
  { 
    return vec3(0,0,0);
  }

  depth += 1;
 
  //Initialise result to be returned
  vec3 result = vec3(0.0f,0.0f,0.0f);

  //============= Triangles =============//
  #pragma region TRIANGLES
  if(current.triangleIndex != -1)
  {
    //============= Smoothness/Mirror =============//
    //If there is any smoothness (for reflectance)
    if(triangles[current.triangleIndex].smoothness == 1.0f && triangles[current.triangleIndex].refractiveIndex == 1.0f)
    {
      //Find reflected ray from indidence ray and normal
      vec3 incidentRay = Vec4ToVec3(current.position) - previous;
      vec3 incidentNorm = normalize(incidentRay);
      vec3 normalVec3 = Vec4ToVec3(triangles[current.triangleIndex].normal);
      vec3 reflectedRay = Reflect(incidentNorm, normalVec3);
      vec4 reflectedRayV4 = Vec3ToHomogenous(reflectedRay);

      //Find the next intersection of the reflected ray
      Intersection nextIntersection;
      bool isIntersection = ClosestIntersection(current.position,reflectedRayV4,triangles,nextIntersection,spheres);

      //If there is an intersection
      if(isIntersection)
      {
        //Add (slight attenuated) directlight from subsequent reflections, taking mirror colour into account
        result += (0.8f * PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres) * triangles[current.triangleIndex].color);
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

    bool castSpecularRay = false;
    //============= Refraction =============//
    if(triangles[current.triangleIndex].refractiveIndex != 1.0f)
    {
      //Random number from 0->1
      float randNum = ((float) rand() / (RAND_MAX));

      if(randNum >= 0.8)
      // if(false)
      {
        castSpecularRay = true;
        // cout << randNum << '\n';
      }
      else
      {
      
      

      //Find refracted ray from indidence ray and normal
      vec3 incidentRay = Vec4ToVec3(current.position) - previous;
      vec4 normalVec4 = triangles[current.triangleIndex].normal;//NormaliseNoHomogenous(current.position - triangles[current.triangleIndex].centre);
      vec3 normalVec3 = Vec4ToVec3(normalVec4);
      vec3 incidentNorm = normalize(incidentRay);
      vec3 refractedRay = Refract(incidentNorm, normalVec3,triangles[current.triangleIndex].refractiveIndex);
      vec4 refractedRayV4 = Vec3ToHomogenous(refractedRay);

      //Find the next intersection of the refracted ray
      Intersection nextIntersection;
      bool isIntersection = ClosestIntersection(current.position,refractedRayV4,triangles,nextIntersection,spheres);

      //If there is an intersection
      if(isIntersection)
      {
        // cout << refractedRayV4 << "\n";
        //Add (slight attenuated) directlight from subsequent reflections, taking mirror colour into account
        result += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres) * triangles[current.triangleIndex].color);
        // cout <<"result" <<result << "\n";
      }
      else//If there is no intersection, don't recurse anymore 
      {
        
      }
      }
    }

    //============= Specularity =============//
    bool castDiffuseRay = false;
    

    if(triangles[current.triangleIndex].smoothness < 1.0f && triangles[current.triangleIndex].smoothness > 0.0f && triangles[current.triangleIndex].refractiveIndex == 1.0f)
    {
      //if 0.6, 60% chance of specular ray being cast, 40% of diffuse
      float smoothness = triangles[current.triangleIndex].smoothness;

      //Random number from 0->1
      float randNum = ((float) rand() / (RAND_MAX));

      if(randNum <= smoothness)//do spec
      {
        //Cast a specular ray, not a diffuse
        // castDiffuseRay = true;
        castSpecularRay = true;

        //============= Random Specular Sample =============//
        //Sample within specular radius
        float specularSampleRadius = 1.0f-smoothness;//0.2f;
        float randXOffset = ((((float) rand() / (RAND_MAX))*specularSampleRadius*2.0f)-specularSampleRadius);
        float randZOffset = ((((float) rand() / (RAND_MAX))*specularSampleRadius*2.0f)-specularSampleRadius);
        vec3 localDirectionOffset = vec3(randXOffset,0,randZOffset);

        //Create local coordinate system of current triangle
        vec3 normal = Vec4ToVec3(triangles[current.triangleIndex].normal);
        vec3 Nt;
        vec3 Nb;
        CreateCoordinateSystem(normal,Nt,Nb);

        //convert to local offset to global coordinates
        vec3 globalDirectionOffset = vec3(localDirectionOffset.x * Nb.x + localDirectionOffset.y * normal.x + localDirectionOffset.z * Nt.x, 
                                          localDirectionOffset.x * Nb.y + localDirectionOffset.y * normal.y + localDirectionOffset.z * Nt.y, 
                                          localDirectionOffset.x * Nb.z + localDirectionOffset.y * normal.z + localDirectionOffset.z * Nt.z);

        //============= Reflected Ray =============//
        //Find reflected ray from indidence ray and normal
        vec3 incidentRay = Vec4ToVec3(current.position) - previous;
        vec3 incidentNorm = normalize(incidentRay);
        vec3 normalVec3 = Vec4ToVec3(triangles[current.triangleIndex].normal);
        vec3 reflectedRay = Reflect(incidentNorm, normalVec3);
        reflectedRay = normalize(reflectedRay);

        //Add offset to reflected ray
        vec3 specularRayVec3 = reflectedRay + globalDirectionOffset;
        vec4 specularRay = Vec3ToHomogenous(specularRayVec3);


        //============= Cast Ray/ Closest Intersection =============//
        //Cast ray
        Intersection nextIntersection;
        bool isIntersection = ClosestIntersection(current.position,specularRay,triangles,nextIntersection,spheres);

        //If there is an intersection
        if(isIntersection)
        {
          // vec3 viewDirection = normalize(previous - Vec4ToVec3(current.position));
          // cout << BRDF << "\n";
          //Add (slight attenuated) directlight from subsequent reflections, taking mirror colour into account
          result += (0.8f * PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres) * triangles[current.triangleIndex].color);
        }
        else//If there is no intersection, don't recurse anymore 
        {
          //(effectively the reflected ray shoots into empty space)
          // result += vec3(0.0f,0.0f,0.0f);
        }
      }
      else//do diffuse
      {
        castDiffuseRay = true;
      }
    }
    //=============End Specularity =============//



    //============= Diffuse =============//
    //If the surface is diffuse (no smoothness)
    if((triangles[current.triangleIndex].smoothness == 0.0f || castDiffuseRay) && triangles[current.triangleIndex].refractiveIndex == 1.0f)
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
      bool isIntersect = ClosestIntersection(current.position,WorldSampleDir,triangles,nextIntersection,spheres);

      //If there is an interesction, recurse at +1 depth
      if(isIntersect)
      {
        // float cosThrefractiveIndexRatio = 2.0f * PI * rand2;
        // indirectLight += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position)) * (cosThrefractiveIndexRatio));

        float cosThrefractiveIndexRatio = rand1;
        indirectLight += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres)*cosThrefractiveIndexRatio*0.5f); //0.5f is the attenuation factor
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
    //============= Direct Lighting =============//
    // if(true)
    // {
    //   // float temp = 2.0f * PI;
    //   // result += ( (DirectLight( current, lightPos, lightColour, triangles ) * triangles[current.triangleIndex].color));
    //   result += ( (AreaLightSample( current, lightPos, lightColour, triangles ) * triangles[current.triangleIndex].color));
    // }
    //In all cases, we sample the light (at maximum importance)

    if(triangles[current.triangleIndex].refractiveIndex == 1.0f || castSpecularRay)
    {
      if(isAreaLight)//arealight
      {
        //if we're casting a specular ray
        if(castSpecularRay)
        {
          //Find reflected ray
          // vec3 incidentRay = Vec4ToVec3(lightPos) - Vec4ToVec3(current.position);
          // vec4 normalVec4 = triangles[current.triangleIndex].normal;//NormaliseNoHomogenous(current.position - triangles[current.triangleIndex].centre);
          // vec3 normalVec3 = Vec4ToVec3(normalVec4);
          // vec3 reflectedRay = normalize(Reflect(incidentRay, normalVec3));

          //Find view direction
          // vec3 viewDirection = normalize(previous - Vec4ToVec3(current.position));

          //Compute direct light
          float BRDF = 1.0f;// pow(abs(glm::dot(viewDirection,reflectedRay)),3); //was 50

          result += ( (AreaLightSample( current, lightPos, lightColour, triangles, spheres ) * triangles[current.triangleIndex].color) * BRDF);
        }
        else //if we are casting a diffuse ray
        {
          
          result += ( (AreaLightSample( current, lightPos, lightColour, triangles, spheres ) * triangles[current.triangleIndex].color));

        }

        // result += ( (AreaLightSample( current, lightPos, lightColour, triangles, spheres ) * triangles[current.triangleIndex].color));
      }
      else//pointlight
      {
        //if we're casting a specular ray
        if(castSpecularRay)
        {
          //Find reflected ray
          // vec3 incidentRay = Vec4ToVec3(lightPos) - Vec4ToVec3(current.position);
          // vec3 incidentNorm = normalize(incidentRay);
          // vec3 normalVec3 = Vec4ToVec3(triangles[current.triangleIndex].normal);
          // vec3 reflectedRay = normalize(Reflect(incidentNorm, normalVec3));

          //Find view direction
          // vec3 viewDirection = normalize(previous - Vec4ToVec3(current.position));

          //Compute direct light
          // cout << (glm::dot(viewDirection,reflectedRay)) << "\n";
          // float BRDF = pow(abs(glm::dot(viewDirection,reflectedRay)),50.0f);
          float BRDF = 1.0f;

          result += ( ((DirectLight( current, lightPos, lightColour, triangles,spheres ) * triangles[current.triangleIndex].color)) * BRDF   );
        }
        else  //if we are casting a diffuse ray
        {
          result += ( (DirectLight( current, lightPos, lightColour, triangles,spheres ) * triangles[current.triangleIndex].color));
        }
      }
    }
    //============= End Direct Lighting =============//
  }
  #pragma endregion TRIANGLES
  //============= End Triangles =============//


  //============= Spheres =============//
  #pragma region SPHERES
  if(current.sphereIndex != -1)
  {
    //============= Smoothness/Mirror =============//
    //If there is any smoothness (for reflectance)
    if(spheres[current.sphereIndex].smoothness == 1.0f && spheres[current.sphereIndex].refractiveIndex == 1.0f)
    {
      //Find reflected ray from indidence ray and normal
      vec3 incidentRay = Vec4ToVec3(current.position) - previous;
      vec4 normalVec4 = NormaliseNoHomogenous(current.position - spheres[current.sphereIndex].centre);
      vec3 normalVec3 = Vec4ToVec3(normalVec4);
      vec3 incidentNorm = normalize(incidentRay);
      vec3 reflectedRay = Reflect(incidentNorm, normalVec3);
      vec4 reflectedRayV4 = Vec3ToHomogenous(reflectedRay);

      //Find the next intersection of the reflected ray
      Intersection nextIntersection;
      bool isIntersection = ClosestIntersection(current.position,reflectedRayV4,triangles,nextIntersection,spheres);

      //If there is an intersection
      if(isIntersection)
      {
        //Add (slight attenuated) directlight from subsequent reflections, taking mirror colour into account
        result += (0.8f * PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres) * spheres[current.sphereIndex].color);
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

    bool castSpecularRay = false;
    //============= Refraction =============//
    if(spheres[current.sphereIndex].refractiveIndex != 1.0f)
    {
      //Random number from 0->1
      // float randNum = ((float) rand() / (RAND_MAX));

      // if(randNum >= 0.8)
      // if()
      // if(false)
      // {
        // castSpecularRay = true;
        // cout << randNum << '\n';
      // }
      // else
      // {
        //Find refracted ray from indidence ray and normal
        vec3 incidentRay = Vec4ToVec3(current.position) - previous;
        vec4 normalVec4 = NormaliseNoHomogenous(current.position - spheres[current.sphereIndex].centre);
        vec3 normalVec3 = Vec4ToVec3(normalVec4);
        vec3 incidentNorm = normalize(incidentRay);

        vec3 refractedRay;// = Refract(incidentNorm, normalVec3,spheres[current.sphereIndex].refractiveIndex);
        float randNum = ((float) rand() / (RAND_MAX));
        if(randNum < 0.9)
        {
          refractedRay = Refract(incidentNorm, normalVec3,spheres[current.sphereIndex].refractiveIndex);
        }
        else
        {
          refractedRay = Reflect(incidentNorm, normalVec3);//,spheres[current.sphereIndex].refractiveIndex);
        }
        
        // vec3 refractedRay;// = Refract(incidentNorm, normalVec3,spheres[current.sphereIndex].refractiveIndex);
        vec4 refractedRayV4 = Vec3ToHomogenous(refractedRay);

        //Find the next intersection of the refracted ray
        Intersection nextIntersection;
        //bool isIntersection = ClosestIntersection(current.position,refractedRayV4,triangles,nextIntersection,spheres);
        bool isIntersection = ClosestIntersection(current.position,refractedRayV4,triangles,nextIntersection,spheres);

        //If there is an intersection
        if(isIntersection)
        {
          result += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres));// * spheres[current.sphereIndex].color);
          // return result;
        }
        else//If there is no intersection, don't recurse anymore 
        {
          
        }
      // }
    }

    //============= End Refraction =============//

    //============= Specularity =============//
    bool castDiffuseRay = false;
    

    if(spheres[current.sphereIndex].smoothness < 1.0f && spheres[current.sphereIndex].smoothness > 0.0f && spheres[current.sphereIndex].refractiveIndex == 1.0f)
    {
      //if 0.6, 60% chance of specular ray being cast, 40% of diffuse
      float smoothness = spheres[current.sphereIndex].smoothness;

      //Random number from 0->1
      float randNum = ((float) rand() / (RAND_MAX));

      if(randNum <= smoothness)//do spec
      {
        //Cast a specular ray, not a diffuse
        // castDiffuseRay = true;
        castSpecularRay = true;

        //============= Random Specular Sample =============//
        //Sample within specular radius
        float specularSampleRadius = 1.0f-smoothness;//0.2f;
        float randXOffset = ((((float) rand() / (RAND_MAX))*specularSampleRadius*2.0f)-specularSampleRadius);
        float randZOffset = ((((float) rand() / (RAND_MAX))*specularSampleRadius*2.0f)-specularSampleRadius);
        vec3 localDirectionOffset = vec3(randXOffset,0,randZOffset);

        //Create local coordinate system of current triangle
        vec4 normalVec4 = NormaliseNoHomogenous(current.position - spheres[current.sphereIndex].centre);
        vec3 normal = Vec4ToVec3(normalVec4);
        vec3 Nt;
        vec3 Nb;
        CreateCoordinateSystem(normal,Nt,Nb);

        //convert to local offset to global coordinates
        vec3 globalDirectionOffset = vec3(localDirectionOffset.x * Nb.x + localDirectionOffset.y * normal.x + localDirectionOffset.z * Nt.x, 
                                          localDirectionOffset.x * Nb.y + localDirectionOffset.y * normal.y + localDirectionOffset.z * Nt.y, 
                                          localDirectionOffset.x * Nb.z + localDirectionOffset.y * normal.z + localDirectionOffset.z * Nt.z);

        //============= Reflected Ray =============//
        //Find reflected ray from indidence ray and normal
        vec3 incidentRay = Vec4ToVec3(current.position) - previous;
        vec3 normalVec3 = normal;
        vec3 incidentNorm = normalize(incidentRay);
        vec3 reflectedRay = Reflect(incidentNorm, normalVec3);
        reflectedRay = normalize(reflectedRay);

        //Add offset to reflected ray
        vec3 specularRayVec3 = reflectedRay + globalDirectionOffset;
        vec4 specularRay = Vec3ToHomogenous(specularRayVec3);


        //============= Cast Ray/ Closest Intersection =============//
        //Cast ray
        Intersection nextIntersection;
        bool isIntersection = ClosestIntersection(current.position,specularRay,triangles,nextIntersection,spheres);

        //If there is an intersection
        if(isIntersection)
        {
          // vec3 viewDirection = normalize(previous - Vec4ToVec3(current.position));
          // float BRDF = pow(abs(glm::dot(viewDirection,reflectedRay)),50.0f);
          // BRDF = 1.0f;

          //Add (slight attenuated) directlight from subsequent reflections, taking mirror colour into account
          result += (0.8f * PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres) * spheres[current.sphereIndex].color);
        }
        else//If there is no intersection, don't recurse anymore 
        {
          //(effectively the reflected ray shoots into empty space)
          // result += vec3(0.0f,0.0f,0.0f);
        }
      }
      else//do diffuse
      {
        castDiffuseRay = true;
      }
    }
    //=============End Specularity =============//



    //============= Diffuse =============//
    //If the surface is diffuse (no smoothness)
    if(spheres[current.sphereIndex].smoothness == 0.0f && spheres[current.sphereIndex].refractiveIndex == 1.0f)// || castDiffuseRay)
    {
      //indirectlight comes from sampling
      vec3 indirectLight;

      //Find local coordinate system at the intersection
      vec4 normalVec4 = NormaliseNoHomogenous(current.position - spheres[current.sphereIndex].centre);
      vec3 normal = Vec4ToVec3(normalVec4);
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
      bool isIntersect = ClosestIntersection(current.position,WorldSampleDir,triangles,nextIntersection,spheres);

      //If there is an interesction, recurse at +1 depth
      if(isIntersect)
      {
        // float cosThrefractiveIndexRatio = 2.0f * PI * rand2;
        // indirectLight += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position)) * (cosThrefractiveIndexRatio));

        float cosThrefractiveIndexRatio = rand1;
        indirectLight += (PathTracer(nextIntersection,lightPos,lightColour,triangles,depth,Vec4ToVec3(current.position),true,isAreaLight,spheres)*cosThrefractiveIndexRatio*0.5f); //0.5f is the attenuation factor
      }
      else//If there is no intersection, do nothing
      {
        //do nothing
      }

      //Div. indirect light by pdf
      indirectLight = (indirectLight/PDF);


      // Backtrace light from the (non-reflective) surface
      result += (indirectLight) * spheres[current.sphereIndex].color;
    }
    //============= End Diffuse =============//



    //============= Emission =============//
    result += spheres[current.sphereIndex].emissive;
    //============= End Emission =============//

    //Direct lighting
    //============= Direct Lighting =============//
    if(spheres[current.sphereIndex].refractiveIndex == 1.0f)// || castSpecularRay)
    {    
      //In all cases, we sample the light (at maximum importance)
      if(isAreaLight)//arealight
      {
        //if we're casting a specular ray
        if(castSpecularRay)
        {
          //Find reflected ray
          vec3 incidentRay = Vec4ToVec3(lightPos) - Vec4ToVec3(current.position);
          vec3 incidentNorm = normalize(incidentRay);
          vec4 normalVec4 = NormaliseNoHomogenous(current.position - spheres[current.sphereIndex].centre);
          vec3 normalVec3 = Vec4ToVec3(normalVec4);
          vec3 reflectedRay = normalize(Reflect(incidentNorm, normalVec3));

          //Find view direction
          vec3 viewDirection = normalize(previous - Vec4ToVec3(current.position));

          //Compute direct light
          float BRDF = pow(abs(glm::dot(viewDirection,reflectedRay)),50.0f);

          result += ( (AreaLightSample( current, lightPos, lightColour, triangles, spheres ) * spheres[current.sphereIndex].color) * BRDF);
        }
        else //if we are casting a diffuse ray
        {
          
          result += ( (AreaLightSample( current, lightPos, lightColour, triangles, spheres ) * spheres[current.sphereIndex].color));

        }

        
      }
      else//pointlight
      {
        //if we're casting a specular ray
        if(castSpecularRay)
        {
          //Find reflected ray
          vec3 incidentRay = Vec4ToVec3(lightPos) - Vec4ToVec3(current.position);
          vec3 incidentNorm = normalize(incidentRay);
          vec4 normalVec4 = NormaliseNoHomogenous(current.position - spheres[current.sphereIndex].centre);
          vec3 normalVec3 = Vec4ToVec3(normalVec4);
          vec3 reflectedRay = normalize(Reflect(incidentNorm, normalVec3));

          //Find view direction
          vec3 viewDirection = normalize(previous - Vec4ToVec3(current.position));

          //Compute direct light
          float BRDF = pow(abs(glm::dot(viewDirection,reflectedRay)),50.0f);

          result += ( ((DirectLight( current, lightPos, lightColour, triangles,spheres ) * spheres[current.sphereIndex].color)) * BRDF   );
        }
        else  //if we are casting a diffuse ray
        {
          
          result += ( (DirectLight( current, lightPos, lightColour, triangles,spheres ) * spheres[current.sphereIndex].color));

        }
      }
    }
    //============= End Direct Lighting =============//
  }


  #pragma endregion SPHERES


  //Return result
  return result;
}

vec3 AreaLightSample( Intersection& intersection, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles,const vector<Sphere>& spheres )
{
  //============= Find Hits to Area Light =============//
  int hits = 0;

  //sample SPHERE_LIGHT_SAMPLES points on disc
  float randX;
  float randY;
  float randZ;

  float threfractiveIndexRatio = 0;
  float phi = 0;
  for(int i = 0; i < SPHERE_LIGHT_SAMPLES; i++)
  {
    //Spherical sampling
    threfractiveIndexRatio = ((float) rand() / (RAND_MAX)) * 2.0f * PI;
    phi = acos((((float) rand() / (RAND_MAX)) * 2.0f) - 1.0f);

    randX = cosf(threfractiveIndexRatio)*sinf(phi);
    randY = sinf(threfractiveIndexRatio)*sinf(phi);
    randZ = cosf(phi);

    randX = (randX*SPHERE_LIGHT_RADIUS*2.0f) - SPHERE_LIGHT_RADIUS;
    randY = (randY*SPHERE_LIGHT_RADIUS*2.0f) - SPHERE_LIGHT_RADIUS;
    randZ = (randZ*SPHERE_LIGHT_RADIUS*2.0f) - SPHERE_LIGHT_RADIUS;

    //Find light/ intersections
    vec4 samplePoint = lightPos + vec4(randX,randY,randZ,0);


    vec4 lightDir = samplePoint - intersection.position;
    vec4 lightDirNormalised = NormaliseNoHomogenous(samplePoint - intersection.position);

    float lightDist = glm::length(lightDir);

    Intersection lightIntersection;
    bool intersectionFound = ClosestIntersection(intersection.position,lightDirNormalised,triangles,lightIntersection,spheres);
    
    //If no intersection found, then we hit the light
    if(!intersectionFound)
    {
      hits+=1;
    }
    //If intersection found, but intersection is in or beyond light, then we hit the light
    else if(glm::length(lightIntersection.position - intersection.position) >= lightDist)
    {
      hits +=1;
    }
  }//end random sample hit counting

  //============= Calcuate power =============//
  float proportionHits = (float)hits / (float)SPHERE_LIGHT_SAMPLES;

  // Light colour is P
  vec3 P = lightColour * proportionHits;
  vec4 nNorm;
  if(intersection.triangleIndex != -1)
  {
    // Get normal to triangle
    nNorm = NormaliseNoHomogenous(triangles[intersection.triangleIndex].normal);
  }
  else if(intersection.sphereIndex != -1)
  {
    nNorm = NormaliseNoHomogenous(intersection.position - spheres[intersection.sphereIndex].centre);
  }

  // r is vector from intersection point to light source
  vec4 r  = lightPos - intersection.position;
  vec4 rNorm = NormaliseNoHomogenous(r);

  //Compute power per area B
  float A = 4.0f * M_PI * ( pow(glm::length(r),2.0f) );

  vec3 B = vec3(P.x/A,P.y/A,P.z/A);
 
  //Compute power per real surface D
  vec3 D = B * max(glm::dot (rNorm,nNorm) , 0.0f);

  return D;
}


vec3 Reflect(const vec3 &incident, const vec3 &normal)
{
  return (incident - ((2.0f * glm::dot(incident,normal)) * normal));
}

vec3 Refract(const vec3 &incidentNormRay, const vec3 &normal, const float &IoR)
{
    float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(incidentNormRay, normal));
    float prevMediumIoR = 1.0f;//index of refraction of previous medium
    float nextMediumIoR = IoR;//index of refraction of medium we are entereing
    
    vec3 n = normal;
    //Check if hitting from outside or inside
    if (cosi < 0.0f)//cosine of the the incident
    {
      cosi = -cosi;
    } 
    else 
    { 
      swap(prevMediumIoR, nextMediumIoR);//swap if we're leaving the refractive geometry
      n = -normal;//reverse normal
      // cout << "reversenormal\n";
    }


    //Check for TotalInternalReflection
    float refractiveIndexRatio = prevMediumIoR / nextMediumIoR;
    float k = 1.0f - (refractiveIndexRatio * refractiveIndexRatio) * (1.0f - (cosi * cosi));
    if(k < 0.0f)
    {
      // total internal reflection. There is no refraction in this case
      return Reflect(incidentNormRay,normal);//this might be wrong
      cout << "TotalInternalReflection\n";
    }
    else//if no T.I.R.
    {
      return (refractiveIndexRatio * incidentNormRay + (refractiveIndexRatio * cosi - sqrtf(k)) * n);
    }
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
    // cos(threfractiveIndexRatio) = r1 = y
    // cos^2(threfractiveIndexRatio) + sin^2(threfractiveIndexRatio) = 1 -> sin(threfractiveIndexRatio) = srtf(1 - cos^2(threfractiveIndexRatio))
    float sinThrefractiveIndexRatio = sqrtf(1.0f - rand1 * rand1);
    float phi = 2 * M_PI * rand2; 
    float x = sinThrefractiveIndexRatio * cosf(phi); 
    float z = sinThrefractiveIndexRatio * sinf(phi); 
    return vec3(x, rand1, z);


  // const float r = sqrtf(rand1);
  // const float threfractiveIndexRatio = 2.0f * PI * rand2;
  //   const float x = r * cosf(threfractiveIndexRatio);
  //   const float y = r * sinf(threfractiveIndexRatio);
 
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
vec4 NormaliseNoHomogenous(vec4 vector4)
{
  vec3 tempVec3 = vec3(vector4.x,vector4.y,vector4.z);
  tempVec3 = normalize(tempVec3);
  return vec4(tempVec3.x,tempVec3.y,tempVec3.z,1.0f);
}

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

