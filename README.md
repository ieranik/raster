This repository contains my implementation of a software emulator of the classic raster based graphics pipeline. The code takes as input a OpenGL-like code snippet from the `scene.txt` file which contains the camera position and orientation, camera parameters, and draw triangle command with the analogs of translate, rotate, scale, push matrix, and pop matrix commands. The code applies the basic transformations, camera view transformation, and projective transformation. The outputs are stored in the `outputs` folder.
   
# Scene Description
The scene description is provided in the `scene.txt` file. Lines 1-3 of `scene.txt` state the parameters of the `gluLookAt` function, i.e., eye position (eyeX, eyeY, and eyeZ in Line 1), look position (lookX, lookY, and lookZ in Line 2), and up direction (upX, upY, and upZ in Line 3). Line 4 provides the `gluPerspective` parameters, i.e., field of view along Y axis (fovY), aspect ratio indicating the ratio between the field of view along X and the field of view along Y axis (aspectRatio), near distance (near), and far distance (far). The rest of `scene.txt` contains the display code to generate/draw the model. The display code contains 7 commands as follows:
1.	`triangle` command – this command is followed by three lines specifying the coordinates of the three points of the triangle to be drawn. The points being p1, p2, and p3, 9 double values, i.e., p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, and p3.z indicate the coordinates. 
2.	`translate` command – this command is followed by 3 double values (tx, ty, and tz) in the next line indicating translation amounts along X, Y, and Z axes. This is equivalent to glTranslatef(tx, ty, tz) in OpenGL.
3.	`scale` command – this command is followed by 3 double values (sx, sy, and sz) in the next line indicating scaling factors along X, Y, and Z axes. This is equivalent to glScalef(sx, sy, sz) in OpenGL.
4.	`rotate` command – this command is followed by 4 double values in the next line indicating the rotation angle in degree (angle) and the components of the vector defining the axis of rotation(ax, ay, and az). This is equivalent to glRotatef(angle, ax, ay, az) in OpenGL.
5.	`push` command – This is equivalent to `glPushMatrix` in OpenGL. 
6.	`pop` command – This is equivalent to `glPopMatrix` in OpenGL.
7.	`end` command – This indicates the end of the display code.

# Modeling Transformation 
In the Modeling transformation phase, the display code in `scene.txt` is parsed, the transformed positions of the points that follow each triangle command are determined, and the transformed coordinates of the points are written in `stage1.txt` file. We maintain a stack **S** of transformation matrices which is manipulated according to the commands given in the display code. 

# View Transformation 
In the view transformation phase, the `gluLookAt` parameters in `scene.txt` are used to generate the view transformation matrix **V**, and the points in `stage1.txt` are transformed by **V** and written in `stage2.txt`. To determine **V**, first determine mutually perpendicular unit vectors **l**, **r**, and **u** from the `gluLookAt` parameters, then apply a translation to move the eye/camera to origin, and finally, apply a rotation such that the **l** aligns with the -Z axis, **r** with X axis, and **u** with Y axis. 

# Projection Transformation 
In the projection transformation phase, the `gluPerspective` parameters in `scene.txt` are used to generate the projection transformation matrix **P**, and the points in `stage2.txt` are transformed by **P** and written in `stage3.txt`. To generate **P**, we first compute the field of view along X (fovX) axis and determine **r** and **t** according to [this link](http://www.songho.ca/opengl/gl_projectionmatrix.html).




