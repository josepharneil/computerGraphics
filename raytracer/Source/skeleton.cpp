#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

#define SCREEN_WIDTH 400
#define SCREEN_HEIGHT 400
#define FULLSCREEN_MODE true
#define PI 3.14159265
// #define isAAOn false

/* * * * * * * * * * * * * * * * * * * * * * *
 *              Global Variables
 * * * * * * * * * * * * * * * * * * * * * * */
int t;
bool quit;


/* * * * * * * * * * * * * * * * * * *
* 	          Overrides
* * * * * * * * * * * * * * * * * * */
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


/* * * * * * * * * * * * * * * * * * * * * * *
 *              Structures
 * * * * * * * * * * * * * * * * * * * * * * */
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


/* * * * * * * * * * * * * * * * * * * * * * *
 *              FUNCTION DEFS
 * * * * * * * * * * * * * * * * * * * * * * */
#pragma region FunctionDefs
void Update(vec4& cameraPos, int& yaw, vec4& lightPos, mat4& cameraMatrix, bool& isAAOn);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, 
                           int& yaw, vec4& lightPos, vec3& lightColour, mat4& cameraMatrix,bool& isAAOn);
bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection);
vec3 DirectLight( const Intersection& intersection, vec4& lightPos, 
                  vec3& lightColour, vector<Triangle>& triangles);
mat4 LookAt(const vec3& from, const vec3& to);

#pragma endregion FunctionDefs


/* * * * * * * * * * * * * * * * * * * * * * *
 *                MAIN
 * * * * * * * * * * * * * * * * * * * * * * */
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
  vec4 cameraPos(0.0f,0.0f,-3.0f,1.0f);
  mat4 cameraMatrix;
  int yaw = 0;

  //Create light source
  vec4 lightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec4 originalLightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec3 lightColour = 14.0f * vec3( 1.0f, 1.0f, 1.0f );

  bool isAAOn = true;

  //Update and draw
  while( !quit ) //NoQuitMessageSDL() )
  {
    Update(cameraPos, yaw, originalLightPos, cameraMatrix, isAAOn);

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

    
    Draw(screen, triangles, cameraPos, yaw, lightPos, lightColour, cameraMatrix, isAAOn);




    SDL_Renderframe(screen);
  }

  //Output
  SDL_SaveImage( screen, "screenshot.bmp" );

  //Kill screen
  KillSDL(screen);
  return 0;

}

/* * * * * * * * * * * * * * * * * * * * * * *
 *                DRAW
 * * * * * * * * * * * * * * * * * * * * * * */
/*Place your drawing here*/
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, 
                          int& yaw, vec4& lightPos, vec3& lightColour, mat4& cameraMatrix, bool& isAAOn)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // FOCAL LENGTH
  float focalLength = SCREEN_WIDTH/2;

  //Indirect lighting approximation
  vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );

  //Instantiate closest intersection
  Intersection closestIntersection;

  //For each pixel
  for (int row = 0; row < SCREEN_HEIGHT; row++)
  {
    for (int col = 0; col < SCREEN_WIDTH; col++)
    {
      if(isAAOn)
      {
        /* * * * * * * * * * * * * * * * * * *
        * 	       Anti-Aliasing
        * * * * * * * * * * * * * * * * * * */

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
        int seed = (row * SCREEN_WIDTH) + col;

        srand(seed);
        
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

        vec3 colourTotal;

        //For each subpixel
        for (int i = 0; i < 4; i++)
        {
          //If there is an intersection
          if (isIntersections[i])
          {
            //Get triangle from triangles
            Triangle intersectedTriangle = triangles[closestIntersections[i].triangleIndex];

            //Compute lighting
            vec3 directLight = DirectLight(closestIntersections[i], lightPos, lightColour, triangles);

            //Get colour of triangle
            vec3 colour = (directLight + indirectLight) * intersectedTriangle.color;

            //Running colour total over whole pixel
            colourTotal += colour;
          }
        }
        colourTotal = colourTotal/4.0f;

        //set to colour of that triangle
        PutPixelSDL(screen, col, row, colourTotal);
      }

      else
      {
        //Compute ray direction
        vec4 direction = vec4(
          col-SCREEN_WIDTH/2,
          row-SCREEN_HEIGHT/2,
          focalLength,
          1.0f);

        //Normalise direction of ray
        direction = normalize(direction);

        //Compute ClosestIntersection
        bool intersect = ClosestIntersection(
          vec4(0,0,0,1),
          direction,
          triangles,
          closestIntersection);

        //If an intersection occurs
        if (intersect)
        {
          //Get triangle from triangles
          Triangle intersectedTriangle = triangles[closestIntersection.triangleIndex];

          //Compute lighting
          vec3 directLight = DirectLight(closestIntersection, lightPos, lightColour, triangles);

          //Get colour of triangle
          vec3 colour = (directLight + indirectLight) * intersectedTriangle.color;

          //set to colour of that triangle
          PutPixelSDL(screen, col, row, colour);
        }
      }
    }
  }
}

/* * * * * * * * * * * * * * * * * * * * * * *
 *                UPDATE
 * * * * * * * * * * * * * * * * * * * * * * */
/*Place updates of parameters here*/
void Update(vec4& cameraPos, int& yaw, vec4& lightPos, mat4& cameraMatrix, bool& isAAOn)
{
  // static int t = SDL_GetTicks();
  /* Compute frame time */
  // int t2 = SDL_GetTicks();
  // float dt = float(t2-t);
  // t = t2;

  float movementSpeed = 0.1f;

  SDL_Event e;

  // cout << "ASDFsdf " << "\n";
  while(SDL_PollEvent(&e))
  {
    // cout << "asdlgjh \n";

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

  // cameraMatrix = //LookAt(vec3(cameraPos.x,cameraPos.y,cameraPos.z), vec3(0,0,0));
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
}

/* * * * * * * * * * * * * * * * * * * * * * *
 *            Other Functions
 * * * * * * * * * * * * * * * * * * * * * * */

bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection)
{
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


vec3 DirectLight( const Intersection& intersection, vec4& lightPos, 
                        vec3& lightColour, vector<Triangle>& triangles)
{
  // Light colour is P
  vec3 P = lightColour;
  // Get normal to triangle
  vec4 nNorm = normalize(triangles[intersection.triangleIndex].normal);
  // r is vector from intersection point to light source
  vec4 r  = lightPos - intersection.position;
  vec4 rNorm = normalize(r);

  //Compute power per area B
  float A = 4 * M_PI * ( pow(glm::length(r),2) );

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
  bool intersectionFound = ClosestIntersection(intersection.position + 1e-4f * lightDirNorm,lightDirNorm,triangles,lightIntersection);
  
  //If an intersection is found
  //2nd part of conjunction is to deal with the light being inside of the wall 
  //(when the light rays are parallel to the wall, i.e., perpendicular to the normal))
  if (intersectionFound && ( abs(glm::dot( nNorm, (intersection.position - lightPos))) >= 0 ) )
  {
    //If the ray is stopped short by another intersection
    if ( glm::length(lightIntersection.position - intersection.position) < lightDist )
    {
      //Set to black (i.e., a shadow)
      D = vec3(0,0,0);
    }
  }

  return D;
}
