
layout(vertices = 1) out;

uniform sampler2DArray htk;

uniform	mat4 m_pvm;

uniform int width,L;
uniform float windSpeed;
uniform float choppyFactor;


uniform float gridSpacing = 5;
uniform vec2 windowSize;
uniform int pixelsPerEdge = 32;
uniform int gridSize = 512;

uniform int useCulling = 0;
uniform int useAdaTess = 1;

in vec4 posV[];

out vec4 posTC[];

#define ID gl_InvocationID
#define  maxPatchSize 64.0

float height(float u, float v) {

	return 1; // this should be the patch maximum height
}


// Checks if a segment is at least partially inside the frustum
// Need to add a little tolerance in here
bool segmentInFrustum(vec4 p1, vec4 p2) {

	float epsilon = +0.1*p1.w;
	if ((p1.x < -p1.w-epsilon && p2.x < -p2.w-epsilon) || (p1.x > p1.w+epsilon && p2.x > p2.w+epsilon) ||
		(p1.z < -p1.w-epsilon && p2.z < -p2.w-epsilon) || (p1.z > p1.w+epsilon && p2.z > p2.w+epsilon))
		return false;
	else 
		return true;
}


// Measures the screen size of segment p1-p2
float screenSphereSize(vec4 p1, vec4 p2) {

	vec4 viewCenter = (p1+p2) * 0.5;
	vec4 viewUp = viewCenter;
	viewUp.y += distance(p1,p2);
	vec4 p1Proj = viewCenter;
	vec4 p2Proj = viewUp;

	vec4 p1NDC, p2NDC;
	p1NDC = p1Proj/p1Proj.w;
	p2NDC = p2Proj/p2Proj.w;
	
	return( clamp(length((p2NDC.xy - p1NDC.xy) * windowSize * 0.5) / (pixelsPerEdge), 
                    1.0, maxPatchSize));
}


void main() {

	vec2 iLevel;
	vec4 oLevel;

	vec4 posTransV[4];
	vec2 pAux;
	vec4 posTCAux[4];
    


	posTC[ID] = posV[ID];

	
	posTCAux[0] = posTC[0];
	posTCAux[1] = posTC[0] + vec4(0.0, 0.0, gridSpacing, 1.0);
	posTCAux[2] = posTC[0] + vec4(gridSpacing,0.0, 0.0, 1.0);
	posTCAux[3] = posTC[0] + vec4(gridSpacing, 0.0, gridSpacing, 1.0);
	
	float chopiness = choppyFactor;// ( 1+ choppyFactor/(1+(exp(-windSpeed+20)/5)));
	
	vec2 tc = vec2(posTCAux[0].xz) / L  ;
	posTransV[0].xz = posTCAux[0].xz - texture(htk, vec3(tc, LAYER_DX_DZ_SX_SZ)).xy * chopiness;
	
	tc = vec2(posTCAux[1].xz)  / L  ;
	posTransV[1].xz = posTCAux[1].xz - texture(htk, vec3(tc, LAYER_DX_DZ_SX_SZ)).xy * chopiness;
	
	tc = vec2(posTCAux[2].xz) / L  ;
	posTransV[2].xz = posTCAux[2].xz - texture(htk, vec3(tc, LAYER_DX_DZ_SX_SZ)).xy * chopiness;
	
	tc = vec2(posTCAux[3].xz)  / L  ;
	posTransV[3].xz = posTCAux[3].xz - texture(htk, vec3(tc, LAYER_DX_DZ_SX_SZ)).xy * chopiness;	
	
	for (int i = 0; i < 4; ++i ) {
		posTransV[i] = m_pvm * vec4(posTransV[i].x, 0, posTransV[i].z, 1.0);
//		posTCAux[i] = m_pvm * vec4(posTCAux[i].x, 0, posTCAux[i].z, 1.0);
	}

	if (useCulling == 0 ||(		
	            segmentInFrustum(posTransV[0], posTransV[1]) ||
				segmentInFrustum(posTransV[0], posTransV[2]) ||
				segmentInFrustum(posTransV[2], posTransV[3]) ||
				segmentInFrustum(posTransV[3], posTransV[1]))) {
					
		if (useAdaTess == 1) {
					
		// Screen size based LOD

			oLevel = vec4(screenSphereSize(posTransV[0], posTransV[1]),
						screenSphereSize(posTransV[0], posTransV[2]),
						screenSphereSize(posTransV[2], posTransV[3]),
						screenSphereSize(posTransV[3], posTransV[1]));
			iLevel = vec2(max(oLevel[1] , oLevel[3]) , max(oLevel[0] , oLevel[2]) );
/*	
			oLevel = vec4(screenSphereSize(posTCAux[0], posTCAux[1]),
						screenSphereSize(posTCAux[0], posTCAux[2]),
						screenSphereSize(posTCAux[2], posTCAux[3]),
						screenSphereSize(posTCAux[3], posTCAux[1]));
			iLevel = vec2(max(oLevel[1] , oLevel[3]) , max(oLevel[0] , oLevel[2]) );
*/				
		}
		else {
			oLevel = vec4(maxPatchSize);
			iLevel = vec2(maxPatchSize);
		
		}
	}
	else if (useCulling == 1) {
		oLevel = vec4(0);
		iLevel = vec2(0);
		
	}
	else {
		oLevel = vec4(maxPatchSize);
		iLevel = vec2(maxPatchSize);
		
	}

	gl_TessLevelOuter[0] = oLevel[0];
	gl_TessLevelOuter[1] = oLevel[1];
	gl_TessLevelOuter[2] = oLevel[2];
	gl_TessLevelOuter[3] = oLevel[3];
	gl_TessLevelInner[0] = iLevel[0];
	gl_TessLevelInner[1] = iLevel[1];
}