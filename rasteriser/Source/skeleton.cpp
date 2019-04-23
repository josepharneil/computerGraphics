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
#define NEAR_CLIP 0.5f
#define FAR_CLIP 5.0f
#define ANGLE_OF_VIEW 3.14159265/2  //field of view is 90 deg as long as focal length is half screen dimension 

//============= Global Variables =============//
bool quit;
bool isWireframe;

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

vector<Triangle> Clip(Triangle& triangle);
vector<Triangle> Triangulate(vector<vec4> vertices, const vec3 color);
float DotNoHomogenous(const vec4 A, const vec4 B);
vector<vec4> ClipToPlane(vector<vec4>& inputVertices, vec4 planePoint, vec4 planeNormal);


#pragma endregion FunctionDefinitions
//============= END Function Definitions =============//

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
//============= END Overrides =============//

//============= Main =============//
int main( int argc, char* argv[] )
{
  //Initially, do not quit
  quit = false;

  // Vertex A = {.position = vec4(0,0,0,1)};
  // Vertex B = {.position = vec4(1,1,1,1)};
  // Vertex C = {.position = vec4(2,2,2,1)};
  // Vertex D = {.position = vec4(3,3,3,1)};
  // Vertex E = {.position = vec4(4,4,4,1)};
  
  // vector<Vertex> myVertices;
  // myVertices.push_back(A);
  // myVertices.push_back(B);
  // myVertices.push_back(C);
  // myVertices.push_back(D);
  // myVertices.push_back(E);
  // Triangulate(myVertices, vec3(1,1,1));

  isWireframe = true;
  
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

      // w divide
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

  //vector to hold all clipped triangles
  vector<Triangle> clippedTriangles;
  

  //Apply clipping to each triangle
  for( uint32_t i=0; i<triangles.size(); ++i )
  {
    //temp holds the resulting triangles from clipping a single triangle
     vector<Triangle> temp = Clip(triangles[i]);

     for( uint32_t j=0; j<temp.size(); ++j )
     {
       //add each triangle from the result vector to clippedTriangles
       clippedTriangles.push_back(temp[j]);
     }
  }
      
  //For each clipped triangle
  for( uint32_t i=0; i<clippedTriangles.size(); ++i )
  {
    //Construct vector of triangle vertices
    vector<Vertex> vertices(3);
    vertices[0].position = clippedTriangles[i].v0;
    vertices[1].position = clippedTriangles[i].v1;
    vertices[2].position = clippedTriangles[i].v2;

    vec4 currentNormal = clippedTriangles[i].normal;
    vec3 currentReflectance = clippedTriangles[i].color;

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

    if ( e.key.keysym.scancode == SDL_SCANCODE_N )
    {
      isWireframe = !isWireframe;
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
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec4& cameraPos, float& focalLength, 
                 float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& lightPos, vec3& lightPower, 
                 vec3& indirectLightPowerPerArea, vec4& currentNormal, vec3& currentReflectance )
{

  //Find number of vertices of polygon (3 for triangle)
  int V = vertices.size();

  //Initialise 3 long vector of projected vertices (pixel locations)
  vector<Pixel> vertexPixels( V );

  //For each vertex
  for( int i=0; i<V; ++i )
  {
    //Compute projection
    PerspectiveProject( vertices[i], vertexPixels[i], cameraPos, focalLength);
  }
  
  //Initialise vectors to store left-most and right-most positions of each row of the projected triangle
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;

  //Calculates leftPixels and rightPixels
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels, screen );

  //Draws the rows
  DrawPolygonRows( screen, leftPixels, rightPixels, depthBuffer, currentNormal, currentReflectance, lightPos, lightPower, indirectLightPowerPerArea );
}


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
  if(isWireframe)
  {  
    for(int i = 0; i < leftPixels.size(); i++)
    {
      PutPixelSDL(screen, leftPixels[i].x, leftPixels[i].y, vec3(1,1,1));
    }
    for(int i = 0; i < rightPixels.size(); i++)
    {
      PutPixelSDL(screen, rightPixels[i].x, rightPixels[i].y, vec3(1,1,1));
    }
  }
  else
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

//clip triangle, and return a vector of the resulting triangles
vector<Triangle> Clip(Triangle& triangle)
{

  //until triangulation, we work with vec4s instead of Vertex structs, for ease

  //vertices is the working list of vertices. We first add the vertices from the input triangle. Then we clip to each plane of the view frustum in turn. Then we triangulate
  vector<vec4> vertices;
  
  //add vertices from input triangle
  vertices.push_back(triangle.v0);
  vertices.push_back(triangle.v1);
  vertices.push_back(triangle.v2);

  //clip to near plane
  

  //clip to far plane
  vertices = ClipToPlane(vertices,vec4(0,0,FAR_CLIP ,1),vec4(0,0,-1,1));

  vertices = ClipToPlane(vertices,vec4(0,0,NEAR_CLIP,1),vec4(0,0,1,1));

  //clip to top
  // vertices = ClipToPlane(vertices,vec4(0,0,FAR_CLIP,1),vec4(0,0,-1,1));

  //clip to left
  // vertices = ClipToPlane(vertices,vec4(0,0,FAR_CLIP,1),vec4(0,0,-1,1));

  //clip to bottom

  //clip to right

  //... (do for all 6 planes of view frustum)
  
  //triangulate the polygon (might already be a triangle)
  vector<Triangle> result = Triangulate(vertices,triangle.color);
  return result;
}

//use fan triangulation to triangluate the given polygon (specified by ordered list of vertices), return a vector of triangles.
//color is only passed so that we can create triangles
vector<Triangle> Triangulate(vector<vec4> vertices, const vec3 color)
{
  vector<Triangle> result;
  
  //if the passed in polygon is a triagle, just place it in the list and return immediately.
  if(vertices.size() == 3)
  {
    result.push_back(Triangle(vertices[0],vertices[1],vertices[2],color));
    return result;
  }
  
  //number of triangles resulting from triangulation
  int numTriangles = ceil(((float)vertices.size()) /2.0f);

  // cout << vertices.size() << "num tri: " << numTriangles << "\n";
  for(int i = 0; i<numTriangles; i++ )
  {
    result.push_back(Triangle(vertices[0],vertices[i+1],vertices[i+2],color));
    // cout << i << "\n";
    // cout << result[i].v0 << "\n";
    // cout << result[i].v1 << "\n";
    // cout << result[i].v2 << "\n";
  }
  
  return result;
}

float DotNoHomogenous(const vec4 A, const vec4 B)
{
  vec3 tempA = vec3(A.x,A.y,A.z);
  vec3 tempB = vec3(B.x,B.y,B.z);
  return glm::dot(tempA, tempB);
}

//clips a convex polygon to the plane specified by the point planePoint and the normal planeNormal, returns a vector<vec4> specifying the resultant polygon
//note that the order of vertices is important to specifying the 
vector<vec4> ClipToPlane(vector<vec4>& inputVertices, vec4 planePoint, vec4 planeNormal)
{
  vector<vec4> result;

  int n = inputVertices.size();

  float pdot = 0;
  float idot = DotNoHomogenous(planeNormal,(inputVertices[0]-planePoint));
  for(int i = 0; i<n;i++)
  {
    float dot = DotNoHomogenous(planeNormal, (inputVertices[i] - planePoint));
    if(dot * pdot < 0)
    {
      float t = pdot/(pdot - dot);
      vec4 I = inputVertices[i-1] + t * (inputVertices[i] - inputVertices[i-1]);
      result.push_back(I);
    }

    if(dot > 0)
    {
      result.push_back(inputVertices[i]);
    }

    pdot = dot;
  }

  if (pdot * idot < 0)
  {
    float t = pdot/(pdot - idot);
    vec4 I = inputVertices[n-1] + t*(inputVertices[0] - inputVertices[n-1]);
    result.push_back(I);
  }
  
  return result;
  
  // return inputVertices;
}

