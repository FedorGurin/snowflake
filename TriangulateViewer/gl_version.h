#ifndef _GL_VERSION_PROFILE_H_
#define _GL_VERSION_PROFILE_H_

#include <QGLFormat>
#include <QSurfaceFormat>
#define _CAN_USE_OPENGL_4_

#ifdef _CAN_USE_OPENGL_4_
#include <QOpenGLFunctions_4_3_Compatibility>
typedef QOpenGLFunctions_4_3_Compatibility Current_OpenGL_Version_Profile;
#else
#include <QOpenGLFunctions_3_3_Core>
typedef QOpenGLFunctions_3_3_Core Current_OpenGL_Version_Profile;
#endif

extern Current_OpenGL_Version_Profile *gl;

inline void setVersionProfile(QSurfaceFormat &fmt)
{
#ifdef _CAN_USE_OPENGL_4_
	fmt.setVersion(4,3);
#else
	fmt.setVersion(3,3);
#endif
	fmt.setProfile(QSurfaceFormat::CompatibilityProfile);
}

#endif