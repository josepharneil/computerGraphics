#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

// Used to describe a triangular surface:
class Triangle
{
public:
	glm::vec4 v0;
	glm::vec4 v1;
	glm::vec4 v2;
	glm::vec4 normal;
	glm::vec3 color;
	glm::vec3 emissive;
	float smoothness;

	Triangle( glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, glm::vec3 color, glm::vec3 emissive, float smoothness )
		: v0(v0), v1(v1), v2(v2), color(color), emissive(emissive), smoothness(smoothness)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
	  glm::vec3 e1 = glm::vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
	  glm::vec3 e2 = glm::vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
	  glm::vec3 normal3 = glm::normalize( glm::cross( e2, e1 ) );
	  normal.x = normal3.x;
	  normal.y = normal3.y;
	  normal.z = normal3.z;
	  normal.w = 1.0;
	}
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles )
{
	using glm::vec3;
	using glm::vec4;

	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );
	vec3 pureWhite(  1.0f, 1.0f, 1.0f );
	vec3 black(  0.00f, 0.00f, 0.00f );

	triangles.clear();
	triangles.reserve( 5*2*3 );

	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	vec4 A(L,0,0,1);
	vec4 B(0,0,0,1);
	vec4 C(L,0,L,1);
	vec4 D(0,0,L,1);

	vec4 E(L,L,0,1);//back-top-left
	vec4 F(0,L,0,1);//front-top-left
	vec4 G(L,L,L,1);//back-top-right
	vec4 H(0,L,L,1);//front-top-right



	float areaLightLength = 60.0f;
	//x,y,z,w
	vec4 areaE(L-areaLightLength, L-0.5f, areaLightLength  , 1);
	vec4 areaF(areaLightLength  , L-0.5f, areaLightLength  , 1);
	vec4 areaG(L-areaLightLength, L-0.5f, L-areaLightLength, 1);
	vec4 areaH(areaLightLength  , L-0.5f, L-areaLightLength, 1);


	// vec4 tri(areaLightLength,L,0,1);
	// vec4 tri1(0,L,L-areaLightLength,1);

	//Area Light
	// triangles.push_back( Triangle( areaE, areaF, areaG, black, pureWhite,0.0f ) );
	// triangles.push_back( Triangle( areaF, areaH, areaG, black, pureWhite,0.0f ) );


	// Floor:
	// triangles.push_back( Triangle( C, B, A, green ) );
	// triangles.push_back( Triangle( C, D, B, green ) );
	triangles.push_back( Triangle( C, B, A, white,black,0.0f ) );
	triangles.push_back( Triangle( C, D, B, white,black,0.0f ) );

	// Left wall
	// triangles.push_back( Triangle( A, E, C, purple ) );
	// triangles.push_back( Triangle( C, E, G, purple ) );
	triangles.push_back( Triangle( A, E, C, red,black,0.0f ) );
	triangles.push_back( Triangle( C, E, G, red,black,0.0f ) );

	// Right wall
	// triangles.push_back( Triangle( F, B, D, yellow ) );
	// triangles.push_back( Triangle( H, F, D, yellow ) );
	triangles.push_back( Triangle( F, B, D, green,black,0.0f ) );
	triangles.push_back( Triangle( H, F, D, green,black,0.0f ) );

	// Ceiling
	// triangles.push_back( Triangle( E, F, G, cyan ) );
	// triangles.push_back( Triangle( F, H, G, cyan ) );
	triangles.push_back( Triangle( E, F, G, white, black,0.0f ) );
	triangles.push_back( Triangle( F, H, G, white, black,0.0f ) );

	// triangles.push_back( Triangle(F,tri,tri1,cyan,black,0.0f) );



	// Back wall
	triangles.push_back( Triangle( G, D, C, white, black,0.0f ) );
	triangles.push_back( Triangle( G, H, D, white, black,0.0f ) );

	// ---------------------------------------------------------------------------
	// Short block

	A = vec4(290,0,114,1);
	B = vec4(130,0, 65,1);
	C = vec4(240,0,272,1);
	D = vec4( 82,0,225,1);
	       
	E = vec4(290,165,114,1);
	F = vec4(130,165, 65,1);
	G = vec4(240,165,272,1);
	H = vec4( 82,165,225,1);

	// // Front
	// triangles.push_back( Triangle(E,B,A,red) );
	// triangles.push_back( Triangle(E,F,B,red) );

	// // Front
	// triangles.push_back( Triangle(F,D,B,red) );
	// triangles.push_back( Triangle(F,H,D,red) );

	// // BACK
	// triangles.push_back( Triangle(H,C,D,red) );
	// triangles.push_back( Triangle(H,G,C,red) );

	// // LEFT
	// triangles.push_back( Triangle(G,E,C,red) );
	// triangles.push_back( Triangle(E,A,C,red) );

	// // TOP
	// triangles.push_back( Triangle(G,F,E,red) );
	// triangles.push_back( Triangle(G,H,F,red) );
	// Front
	triangles.push_back( Triangle(E,B,A,white,black,0.0f) );
	triangles.push_back( Triangle(E,F,B,white,black,0.0f) );

	// Front
	triangles.push_back( Triangle(F,D,B,white,black,0.0f) );
	triangles.push_back( Triangle(F,H,D,white,black,0.0f) );

	// BACK
	triangles.push_back( Triangle(H,C,D,white,black,0.0f) );
	triangles.push_back( Triangle(H,G,C,white,black,0.0f) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,white,black,0.0f) );
	triangles.push_back( Triangle(E,A,C,white,black,0.0f) );

	// TOP
	triangles.push_back( Triangle(G,F,E,white,black,0.0f) );
	triangles.push_back( Triangle(G,H,F,white,black,0.0f) );
	// ---------------------------------------------------------------------------
	// Tall block

	A = vec4(423,0,247,1);
	B = vec4(265,0,296,1);
	C = vec4(472,0,406,1);
	D = vec4(314,0,456,1);
	       
	E = vec4(423,330,247,1);
	F = vec4(265,330,296,1);
	G = vec4(472,330,406,1);
	H = vec4(314,330,456,1);

	// // Front
	// triangles.push_back( Triangle(E,B,A,blue) );
	// triangles.push_back( Triangle(E,F,B,blue) );

	// // Front
	// triangles.push_back( Triangle(F,D,B,blue) );
	// triangles.push_back( Triangle(F,H,D,blue) );

	// // BACK
	// triangles.push_back( Triangle(H,C,D,blue) );
	// triangles.push_back( Triangle(H,G,C,blue) );

	// // LEFT
	// triangles.push_back( Triangle(G,E,C,blue) );
	// triangles.push_back( Triangle(E,A,C,blue) );

	// // TOP
	// triangles.push_back( Triangle(G,F,E,blue) );
	// triangles.push_back( Triangle(G,H,F,blue) );
	// Front
	triangles.push_back( Triangle(E,B,A,white,black,0.9f) );
	triangles.push_back( Triangle(E,F,B,white,black,0.9f) );

	// Front
	triangles.push_back( Triangle(F,D,B,white,black,0.9f) );
	triangles.push_back( Triangle(F,H,D,white,black,0.9f) );

	// BACK
	triangles.push_back( Triangle(H,C,D,white,black,0.9f) );
	triangles.push_back( Triangle(H,G,C,white,black,0.9f) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,white,black,0.9f) );
	triangles.push_back( Triangle(E,A,C,white,black,0.9f) );

	// TOP
	triangles.push_back( Triangle(G,F,E,white,black,0.9f) );
	triangles.push_back( Triangle(G,H,F,white,black,0.9f) );


	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec4(1,1,1,1);
		triangles[i].v1 -= vec4(1,1,1,1);
		triangles[i].v2 -= vec4(1,1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].v0.w = 1.0;
		triangles[i].v1.w = 1.0;
		triangles[i].v2.w = 1.0;
		
		triangles[i].ComputeNormal();
	}
}

#endif
