#version 430

uniform sampler2D pictureTex;

out vec2 pos;
out vec2 texCoord;

void main()
{
	pos.x = gl_VertexID%2;
	pos.y = gl_VertexID/2;
	texCoord = pos;
	gl_Position = vec4(pos.xy*2-1, 0, 1);
}