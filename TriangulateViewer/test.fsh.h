#version 430

uniform sampler2D pictureTex;

in vec2 texCoord;
out vec4 outputColor;

void main()
{
	vec4 texColor=texture(pictureTex, texCoord);
	outputColor = texColor;
}