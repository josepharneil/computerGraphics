#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;
using glm::vec2;

//============= Defines =============//
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800//256
#define FULLSCREEN_MODE false
#define PI 3.14159265
#define LIGHT_POWER 3.5f
#define NEAR_CLIP 0.001
#define FAR_CLIP 1000
#define FIELD_OF_VIEW 90 

//============= Global Variables =============//
bool quit;

//============= Structures =============//
// struct Triangle
// {
//   vec4 v0;
//   vec4 v1;
//   vec4 v2;
//   vec4 normal;
//   vec3 color;
// };

struct Pixel
{
  int x;
  int y;
  float zinv;
  // vec3 illumination;
  vec4 pos3d;
};

struct Vertex
{
  vec4 position;
  // vec4 normal;
  // vec3 reflectance;
};

//============= Function Definitions =============//
#pragma region FunctionDefinitions
void Update(vec4& cameraPos, int& yaw, mat4& cameraMatrix, vec4& lightPos);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, float& focalLength, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea);
void TestComputePolygonRows(screen* screen);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);

//pixel versions
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result);
void DrawPolygonRows( screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& currentNormal, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea );
void FindLine( Pixel a, Pixel b, vector<Pixel>& lineToDraw);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, screen* screen );
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec4& cameraPos, float& focalLength, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, vec4& currentNormal, vec3& currentReflectance );
void PerspectiveProject( const Vertex& vertex, Pixel& p, vec4& cameraPos, float& focalLength );
void PixelShader(const Pixel& p, screen* s, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& currentNormal, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea );
vec4 NormaliseNoHomogenous(vec4 vector4);
void CalculateCameraMatrix(vec4& camPos, int& yaw, mat4& camMatrix);
void CalculateProjectionMatrix(mat4& projectionMatrix);
#pragma endregion FunctionDefinitions
//============= END Function Definitions =============//

//============= Main =============//
int main( int argc, char* argv[] )
{
  //Initially, do not quit
  quit = false;
  
  // Initialise screen
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // TestComputePolygonRows(screen);

  //Instantiate vector of triangles
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  vector<Triangle> originalTriangles;
  LoadTestModel( originalTriangles );

  //Camera position (fixed)
  vec4 cameraPos( 0, 0, -3.001,1 );
  mat4 cameraMatrix;
  int yaw = 0;
  CalculateCameraMatrix(cameraPos, yaw, cameraMatrix);

  //Light
  vec4 lightPos(0,-0.5,-0.7,1);
  vec4 originalLightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec3 lightPower = LIGHT_POWER*vec3( 1, 1, 1 );
  vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

  //Focal length
  float focalLength = SCREEN_WIDTH/2;

  while( !quit ) //NoQuitMessageSDL() )
    {
      Update(cameraPos, yaw, cameraMatrix, originalLightPos);

      mat4 invCameraMatrix = glm::inverse(cameraMatrix);
      
      for (int t = 0; t < triangles.size(); t++)
      {
        triangles[t].v0 = invCameraMatrix * originalTriangles[t].v0;
        triangles[t].v1 = invCameraMatrix * originalTriangles[t].v1;
        triangles[t].v2 = invCameraMatrix * originalTriangles[t].v2;
        triangles[t].ComputeNormal();
      }

      lightPos = invCameraMatrix * originalLightPos;

      Draw(screen, triangles, cameraPos, focalLength, lightPos, lightPower, indirectLightPowerPerArea);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}
//============= End Main =============//

//============= Draw =============//
/*Place your drawing here*/
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, float& focalLength, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  //Depth buffer
  float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

  for( int y=0; y<SCREEN_HEIGHT; ++y )
  {
    for( int x=0; x<SCREEN_WIDTH; ++x )
    {
      depthBuffer[y][x] = 0;
    }
  }
    
  
  //For each triangle
  for( uint32_t i=0; i<triangles.size(); ++i )
  {
    //Construct vector of triangle vertices
    vector<Vertex> vertices(3);
    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;

    vec4 currentNormal = triangles[i].normal;
    vec3 currentReflectance = triangles[i].color;

    //Draw polygon for each triangle
    DrawPolygon( screen, vertices, cameraPos, focalLength, depthBuffer, 
    lightPos, lightPower, indirectLightPowerPerArea, currentNormal, currentReflectance );

  }
}
//============= END Draw =============//

/*
For all triangles
  DrawPolygon -> compute and draw polygon rows
*/

//============= Update =============//
void Update(vec4& cameraPos, int& yaw, mat4& cameraMatrix, vec4& lightPos)
{
  float movementSpeed = 0.1f;

  SDL_Event e;

  while(SDL_PollEvent(&e))
  {

    if( e.type != SDL_KEYDOWN) {continue;}


    //Move camera with arrow keys
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
      cameraPos.x -= movementSpeed;
    }
    if( e.key.keysym.scancode == SDL_SCANCODE_RIGHT )
    {
      //Move camera right
      cameraPos.x += movementSpeed;
    }
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
  
    CalculateCameraMatrix(cameraPos, yaw, cameraMatrix);
  }
}

/////////
//This draws the polygon from some vertices
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec4& cameraPos, float& focalLength, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, vec4& currentNormal, vec3& currentReflectance )
{

  //Find number of vertices of polygon (3 for triangle)
  int V = vertices.size();

  //initialise the pixel structs
  vector<Pixel> vertexPixels( V );


  //make a copy of vertices because we don't want to alter the construction of the scene.
  vector<Vertex> verticesCopy(3);
  verticesCopy = vertices;


  //CLIPPING START///
  #pragma region Clipping

  //A vertex simply contains a vec4 position. The passed in vertices are in camera space. 
  //We need to enforce the w dimension of each to be 1 (since strictly the homogenous component should be introduced here). 
  //Note that we need to fill three "pixel" structs with the correct information
  for( int i=0; i<V; ++i )
  {
    verticesCopy[i].position.w = 1.0f;
  }

  //Now we need to create the projection 4 x 4 matrix 
  mat4 projectionMatrix;
  CalculateProjectionMatrix(projectionMatrix);

  //Now we need to multiply the vertices by this matrix which will transform the vertices into clip space
  for( int i=0; i<V; ++i )
  {
    verticesCopy[i].position = projectionMatrix * verticesCopy[i].position;
  }
  //Now the vertices in verticesCopy are in clip space and we need to carry out clipping

      ////ACTUAL CLIPPING START////




      /////ACTUAL CLIPPING END/// /


  //Now we need to perform the perspective divide (division by W). Then we are in NDC space, where all three dimensions are in range [-1,1]. --Apparently NDC space is a cartesian space--
  for( int i=0; i<V; ++i )
  {
    verticesCopy[i].position.x = verticesCopy[i].position.x / verticesCopy[i].position.w; 
    verticesCopy[i].position.y = verticesCopy[i].position.y / verticesCopy[i].position.w;
    verticesCopy[i].position.z = verticesCopy[i].position.z / verticesCopy[i].position.w;
    verticesCopy[i].position.w = 1.0f; 
  }

  //set the pos3d attribute and zinv attribute of the pixels:
  for( int i=0; i<V; ++i )
  {
    vertexPixels[i].pos3d = verticesCopy[i].position;
    if(verticesCopy[i].position.z == 0)
    {
      vertexPixels[i].zinv = 0;
    }
    else
    {
      vertexPixels[i].zinv = 1.0f/verticesCopy[i].position.z;
    }
  }

  //Now we need to do a viewport transform to get the point into pixel space (also known as raster space). 
  for( int i=0; i<V; ++i )
  {
    vertexPixels[i].x = (verticesCopy[i].position.x + 1) * 0.5f * (SCREEN_WIDTH  - 1);
    vertexPixels[i].y = (verticesCopy[i].position.y + 1) * 0.5f * (SCREEN_HEIGHT - 1);
  }

  
  #pragma endregion Clipping
  ///CLIPPING END///


  /* OLD PERSPECTIVE
  //Initialise 3 long vector of projected vertices (pixel locations)
  vector<Pixel> vertexPixels( V );

  //For each vertex
  for( int i=0; i<V; ++i )
  {
    //Compute projection
    PerspectiveProject( vertices[i], vertexPixels[i], cameraPos, focalLength);
  }*/

  //CLIPPING BE DONE BY HERE
  

  //Initialise vectors to store left-most and right-most positions of each row of the projected triangle
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;

  //Calculates leftPixels and rightPixels
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels, screen );

  //Draws the rows
  DrawPolygonRows( screen, leftPixels, rightPixels, depthBuffer, currentNormal, currentReflectance, lightPos, lightPower, indirectLightPowerPerArea );
}

/*
//should pass an empty vector of vertices for result
//also need a version for Y and Z
void ClipPolygonComponentX(vector<Vertex>& vertices,float componentFactor, vector<Vertex>& result)
{
  Vertex previousVertex = vertices[vertices.length -1];
  float previousComponent = previousVertex.position.x * componentFactor;
  bool previousInside = previousComponent <= previousVertex.position.w;

  for(int i = 0; i++; i < vertices.length)
  {
    Vertex currentVertex = vertices[i];
    float currentComponent = currentVertex.position.x * componentFactor;
    boolean currentInside = currentComponent <= currentVertex.position.w;

    if(currentInside ^ previousInside) //^ is XOR
    {
      float lerpAmt = (previousVertex.position.w - previousComponent) /
					((previousVertex.position.w - previousComponent) - 
					 (currentVertex.position.w - currentComponent));
    }

    if(currentInside) 
    {
      result.insert(0,currentVertex); //add currentVertex to result 
    }

    previousVertex = currentVertex;
    previousComponent = currentComponent;
    previousInside = currentInside;
  }
}

// componentIndex is 0 for x, 1 for y, 2 for z
boolean ClipPolygonAxis(vector<Vertex> vertices, vector<Vertex> auxilliaryList,int componentIndex)
{

}
*/

//Projects scene point (triangle vertex) onto image plane
void PerspectiveProject( const Vertex& vertex, Pixel& p, vec4& cameraPos, float& focalLength)
{
  //Translates vertex so that camera is at origin, relative to the vertex
  vec4 vertexTransformed = vertex.position; //- cameraPos;

  //vertex[0] is X
  //vertex[1] is Y
  //vertex[2] is Z

  int x = (focalLength * (vertexTransformed[0] / vertexTransformed[2])) + (SCREEN_WIDTH /2);
  int y = (focalLength * (vertexTransformed[1] / vertexTransformed[2])) + (SCREEN_HEIGHT/2);

  //float zinv = 1.0f/glm::distance(cameraPos, vertex);
  float zinv;
  if(vertexTransformed.z == 0)
  {
    zinv = 0;
  }
  else
  {
    zinv = 1.0f/(vertexTransformed.z);
  }

  p.x = x;
  p.y = y;
  p.zinv = zinv;
  p.pos3d = vertex.position;
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, screen* screen )
{
  // 1. Find max and min y-value of the polygon
  //    and compute the number of rows it occupies.
  int yMin;
  int yMax;

  //Initialise two bounds
  yMin = +numeric_limits<int>::max();
  yMax = -numeric_limits<int>::max();

  //For each vertex pixel
  for (int i = 0; i < vertexPixels.size(); i++)
  {
    //Update min and max bounds
    if (vertexPixels[i].y < yMin) { yMin = vertexPixels[i].y;}
    if (vertexPixels[i].y > yMax) { yMax = vertexPixels[i].y; }
  }

  //Number of rows of pixels in the triangle
  int numberOfRows = (yMax - yMin + 1);


  // 2. Resize leftPixels and rightPixels
  //    so that they have an element for each row.
  leftPixels.resize(numberOfRows);
  rightPixels.resize(numberOfRows);

  // 3. Initialize the x-coordinates in leftPixels
  //    to some really large value and the x-coordinates
  //    in rightPixels to some really small value.
  for( int i = 0; i < numberOfRows; ++i )
  {
    leftPixels[i].x  = +numeric_limits<int>::max();
    rightPixels[i].x = -numeric_limits<int>::max();

    leftPixels[i].y  = yMin + i;
    rightPixels[i].y = yMin + i;
  }

  // 4. Loop through all edges of the polygon and use
  //    linear interpolation to find the x-coordinate for
  //    each row it occupies. Update the corresponding
  //    values in rightPixels and leftPixels.
  vector<vector<Pixel> > polygonEdges(3);

  FindLine(vertexPixels[0], vertexPixels[1], polygonEdges[0]);
  FindLine(vertexPixels[1], vertexPixels[2], polygonEdges[1]);
  FindLine(vertexPixels[2], vertexPixels[0], polygonEdges[2]);

  //For each edge
  for ( int i = 0; i < 3; i++)
  {
    vector<Pixel>& edge = polygonEdges[i];
    //For each pixel on edge
    for (int j = 0; j < edge.size(); j++)
    {
      //wireframe test
      // PutPixelSDL(screen, edge[j].x, edge[j].y, vec3(1,1,1));

      Pixel pixel = edge[j];

      int rowIndex = (pixel.y - yMin);

      if(pixel.x < leftPixels[rowIndex].x) 
      { 
        leftPixels[rowIndex].x = pixel.x;
        leftPixels[rowIndex].zinv = pixel.zinv;
        leftPixels[rowIndex].pos3d = pixel.pos3d;
      }
      if(pixel.x > rightPixels[rowIndex].x) 
      { 
        rightPixels[rowIndex].x = pixel.x; 
        rightPixels[rowIndex].zinv = pixel.zinv;
        rightPixels[rowIndex].pos3d = pixel.pos3d;
      }
    }
  }
}

void DrawPolygonRows( screen* screen, 
                      const vector<Pixel>& leftPixels, 
                      const vector<Pixel>& rightPixels, 
                      float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH],
                      vec4& currentNormal, vec3& currentReflectance, 
                      vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea)
{
  vector<Pixel> pixelsBetween(SCREEN_HEIGHT);

  //For each entry in left/right pixels
  for (int i = 0; i < leftPixels.size(); i++)
  {
    int numberOfPixelsBetween = rightPixels[i].x - leftPixels[i].x + 1;

    pixelsBetween.resize(numberOfPixelsBetween);
    InterpolatePixel(leftPixels[i], rightPixels[i], pixelsBetween);

    for(int p = 0; p < pixelsBetween.size(); p++)
    {
      if(pixelsBetween[p].y < SCREEN_HEIGHT && pixelsBetween[p].y >= 0 && pixelsBetween[p].x < SCREEN_WIDTH && pixelsBetween[p].x >=0)
      {
        PixelShader(pixelsBetween[p], screen, depthBuffer, currentNormal, currentReflectance, lightPos, lightPower, indirectLightPowerPerArea );
      } 
    }
  }
}


void TestComputePolygonRows(screen* screen)
{
  vector<Pixel> vertexPixels(3);

  vertexPixels[0].x = 10;
  vertexPixels[0].y = 5;
  vertexPixels[0].zinv = 0.0f;

  vertexPixels[1].x = 5;
  vertexPixels[1].y = 10;
  vertexPixels[1].zinv = 0.0f;

  vertexPixels[1].x = 5;
  vertexPixels[1].y = 10;
  vertexPixels[1].zinv = 0.0f;

  // vertexPixels[2] = {15,15, 0.0f};
  vertexPixels[2].x = 15;
  vertexPixels[2].y = 15;
  vertexPixels[2].zinv = 0.0f;

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels, screen );
  
  for( int row=0; row<leftPixels.size(); ++row )
  {
    cout << "Start: ("
    << leftPixels[row].x << ","
    << leftPixels[row].y << "). "
    << "End: ("
    << rightPixels[row].x << ","
    << rightPixels[row].y << "). " << endl;
  }
}


//Could be sped up with Bresenham's line algorithm
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result)
{
  int N = result.size();

  //Calc steps
  float xStep = (b.x-a.x) / float(max(N-1,1));
  float yStep = (b.y-a.y) / float(max(N-1,1));
  float zStep = (b.zinv-a.zinv) / float(max(N-1,1));
  //The quantity we interpolate is pos3dStep over Z - this accounts for projection
  vec4 pos3dStep = ( (b.pos3d*b.zinv) - (a.pos3d*a.zinv) ) / float(max(N-1,1));

  //Initialise
  float xCurrent = float(a.x);
  float yCurrent = float(a.y);
  float zCurrent = a.zinv;
  vec4 pos3dCurrent( a.pos3d * a.zinv );

  //Interpolate
  for( int i=0; i<N; ++i )
  {
    result[i].x = round(xCurrent);
    result[i].y = round(yCurrent);
    result[i].zinv = zCurrent;
    result[i].pos3d = pos3dCurrent / zCurrent;
    
    xCurrent += xStep;
    yCurrent += yStep;
    zCurrent += zStep;
    pos3dCurrent += pos3dStep;
  }
}



//Finds a line on screen between a and b of some colour
void FindLine( Pixel a, Pixel b, vector<Pixel>& lineToDraw)
{
  // Find difference between two vertices.
  Pixel delta;
  delta.x = glm::abs( a.x - b.x );
  delta.y = glm::abs( a.y - b.y );
  // delta.zinv = glm::abs( a.zinv - b.zinv );//not used
   
  //Find number of pixels needed.
  int numberOfPixels = glm::max( delta.x, delta.y ) + 1;

  //Resize lineToDraw for interpolations
  lineToDraw.resize(numberOfPixels);
  
  //Interpolate between two vertices
  InterpolatePixel( a, b, lineToDraw );
}

void PixelShader(const Pixel& p, screen* screen, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], 
                vec4& currentNormal, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea )
{  
  if(p.zinv > depthBuffer[p.y][p.x])
  {
    vec4 r = lightPos - p.pos3d;
    vec4 n = currentNormal;

    vec4 rNorm = NormaliseNoHomogenous(r);
    vec4 nNorm = NormaliseNoHomogenous(n);
    
    float A = 4 * M_PI * ( pow(glm::length(r),2) );
    vec3 B = vec3(lightPower.x/A,lightPower.y/A,lightPower.z/A);

    vec3 D = B * max(glm::dot (rNorm,nNorm), 0.0f);

    vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea);

    PutPixelSDL(screen, p.x, p.y, illumination);

    depthBuffer[p.y][p.x] = p.zinv;
  }
}

//Helpful Functions
vec4 NormaliseNoHomogenous(vec4 vector4)
{
  vec3 tempVec3 = vec3(vector4.x,vector4.y,vector4.z);
  tempVec3 = normalize(tempVec3);
  return vec4(tempVec3.x,tempVec3.y,tempVec3.z,1.0f);
}

void CalculateCameraMatrix(vec4& camPos, int& yaw, mat4& camMatrix)
{
  camMatrix[0][0] = cos( yaw * PI / 180 );
  camMatrix[0][1] = 0;
  camMatrix[0][2] = sin( yaw * PI / 180 );
  camMatrix[0][3] = 0;

  camMatrix[1][0] = 0;
  camMatrix[1][1] = 1;
  camMatrix[1][2] = 0;
  camMatrix[1][3] = 0;

  camMatrix[2][0] = -sin( yaw * PI / 180 );
  camMatrix[2][1] = 0;
  camMatrix[2][2] = cos( yaw * PI / 180 );
  camMatrix[2][3] = 0;

  camMatrix[3][0] = camPos.x;
  camMatrix[3][1] = camPos.y;
  camMatrix[3][2] = camPos.z;
  camMatrix[3][3] = 1;
}

void CalculateProjectionMatrix(mat4& projectionMatrix)
{
  float S = 1.0f/tan((FIELD_OF_VIEW/2.0f)*(PI / 180.0f));

  projectionMatrix[0][0] = S;
  projectionMatrix[0][1] = 0.0f;
  projectionMatrix[0][2] = 0.0f;
  projectionMatrix[0][3] = 0.0f;

  projectionMatrix[1][0] = 0.0f;
  projectionMatrix[1][1] = S;
  projectionMatrix[1][2] = 0.0f;
  projectionMatrix[1][3] = 0.0f;

  projectionMatrix[2][0] = 0.0f;
  projectionMatrix[2][1] = 0.0f;
  projectionMatrix[2][2] = -FIELD_OF_VIEW/(FIELD_OF_VIEW - NEAR_CLIP);
  projectionMatrix[2][3] = -1.0f;

  projectionMatrix[3][0] = 0.0f;
  projectionMatrix[3][1] = 0.0f;
  projectionMatrix[3][2] = -(FIELD_OF_VIEW * NEAR_CLIP)/(FIELD_OF_VIEW - NEAR_CLIP);
  projectionMatrix[3][3] = 0.0f;
}
