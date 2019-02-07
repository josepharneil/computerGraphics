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

#define SCREEN_WIDTH 64   //320
#define SCREEN_HEIGHT 64  //256
#define FULLSCREEN_MODE true

// struct Triangle
// {
//   vec4 v0;
//   vec4 v1;
//   vec4 v2;
//   vec4 normal;
//   vec3 color;
// };
/* GLOBAL VARIABLES                                                      */
int t;
bool quit;

struct Intersection
{
    vec4 position;
    float distance;
    int triangleIndex;
};

/* * * * * * * * * * * * * * * * * * * * * * *
 *              FUNCTION DEFS
 * * * * * * * * * * * * * * * * * * * * * * */

void Update(vec4& cameraPos, float& yaw);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, 
                           float& yaw, vec4& lightPos, vec3& lightColour);
bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection);
vec3 DirectLight( const Intersection& i, vec4& lightPos, vec3& lightColour, vector<Triangle>& triangles);



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

  //Camera control
  vec4 cameraPos(0.0f,0.0f,-3.0f,1.0f);
  float yaw = 0.0f;

  //Create light source
  vec4 lightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec3 lightColour = 14.0f * vec3( 1.0f, 1.0f, 1.0f );

  //Update and draw
  while( 1 ) //NoQuitMessageSDL() )
    {
      Draw(screen, triangles, cameraPos, yaw, lightPos, lightColour);
      Update(cameraPos, yaw);
      SDL_Renderframe(screen);

      if(quit){break;}

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
                          float& yaw, vec4& lightPos, vec3& lightColour)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // Create rotation matrix
  mat4 rotMat( vec4(1.0f,0.0f,0.0f,yaw), 
               vec4(0.0f,1.0f,0.0f,0.0f),
               vec4(-yaw,0.0f,1.0f,0.0f),
               vec4(0.0f,0.0f,0.0f,1.0f) );

   vec4 right(   rotMat[0][0], rotMat[0][1], rotMat[0][2], 1 );
   vec4 down(    rotMat[1][0], rotMat[1][1], rotMat[1][2], 1 );
   vec4 forward( rotMat[2][0], rotMat[2][1], rotMat[2][2], 1 );


  float focalLength = 20.0f;//250.0f;//bigger zooms

  Intersection closestIntersection;

  //For each pixel
  for (int row = 0; row < SCREEN_HEIGHT; row++)
  {
    for (int col = 0; col < SCREEN_WIDTH; col++)
    {
      //Compute ray direction
      vec4 direction = vec4(
        col-SCREEN_WIDTH/2,
        row-SCREEN_HEIGHT/2,
        focalLength,
        1.0f);

      direction = normalize(direction);

      //Compute ClosestIntersection
      bool intersect = ClosestIntersection(
        cameraPos,
        direction,
        triangles,
        closestIntersection);

      //If an intersection occurs
      if (intersect)
      {
        //Get triangle from triangles
        Triangle intersectedTriangle = triangles[closestIntersection.triangleIndex];

        //Get colour of triangle
        vec3 colour = intersectedTriangle.color;

        //set to colour of that triangle
        PutPixelSDL(screen, col, row, colour);
      }
      //ELSE set to black
    }
  }
}

/* * * * * * * * * * * * * * * * * * * * * * *
 *                UPDATE
 * * * * * * * * * * * * * * * * * * * * * * */
/*Place updates of parameters here*/
void Update(vec4& cameraPos, float& yaw)
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  float movementSpeed = 0.1f;

  SDL_Event e;

  // cout << "ASDFsdf " << "\n";
  while(SDL_PollEvent(&e))
  {
    // cout << "asdlgjh \n";

    if( e.type != SDL_KEYDOWN) {continue;}

    if( e.key.keysym.scancode == SDL_SCANCODE_UP )
    {
      cameraPos.z += movementSpeed;

    }
    if( e.key.keysym.scancode == SDL_SCANCODE_DOWN )
    {
      //Move camera backward
      cameraPos.z -= movementSpeed;

    }
    if( e.key.keysym.scancode == SDL_SCANCODE_LEFT )
    {
      //Move camera left
      // cameraPos.x -= movementSpeed;

      //Rotate camera anticlockwise around y axis
      yaw -= 10.0f;// % 360;
    }
    if( e.key.keysym.scancode == SDL_SCANCODE_RIGHT )
    {
      //Move camera right
      // cameraPos.x += movementSpeed;


      //Rotate camera clockwise around y axis
      yaw += 10.0f;// % 360;
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

}

/* * * * * * * * * * * * * * * * * * * * * * *
 *               Other Functions
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

    // cout << (decltype(triangle.v0.x)) << " ";

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
    
    // cout << -dir3.x << " ";

    //x = (t u v)^T
    //COME BACK HERE todo cramar's rule instead of inbuilt inverse
    vec3 x = glm::inverse( A ) * b;

    float t = x.x;//if neg, can reject immediately
    float u = x.y;
    float v = x.z;

    // cout << t << " ";
    // if (u < 1 && u > 0) {cout << u << " ";}

    //Check if is intersection is within triangle
    if ( (u >= 0.0f) && (v >= 0.0f) && (u + v <= 1.0f) && (t > 0.0f) )
    {
      result = true;
      //Compute "intersection" structure
      if (t < minDist)
      {
        //We set w = 1 -- not sure of this
        closestIntersection.position = vec4(t*dir.x, t*dir.y, t*dir.z,1);
        closestIntersection.distance = t;
        closestIntersection.triangleIndex = i;

        minDist = t;
      }
    }
  }
  return result;
}



vec3 DirectLight( const Intersection& intersection, vec4& lightPos, vec3& lightColour, vector<Triangle>& triangles)
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
}