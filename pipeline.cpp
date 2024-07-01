#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <stack>
#include <string.h>
#include <vector>

#define pi (2*acos(0.0))
#define eps 0.0000001

using namespace std;

class point
{
public:
    double x,y,z;
    point(double xx=0.0, double yy=0.0, double zz=0.0)
    {
        x=xx;
        y=yy;
        z=zz;
    }
};


class hPoint
{
public:
    double p[4];
    hPoint(double x=0.0, double y=0.0, double z=0.0)
    {
        p[0]=x;
        p[1]=y;
        p[2]=z;
        p[3]=1.0;
    }
    void normalize()
    {
        if(fabs(p[3]-1.0)>eps)
        {
            p[0]/=p[3];
            p[1]/=p[3];
            p[2]/=p[3];
            p[3]=1.0;
        }
    }
};


class vec
{
public:
	double x,y,z;
	vec(){x=0.0;y=0.0;z=0.0;}
	vec(double xx,double yy,double zz){x=xx;y=yy;z=zz;}
	void setq(double xx,double yy,double zz){x=xx;y=yy;z=zz;}
	void norm()
	{
		double len=sqrt(x*x+y*y+z*z);
		x=x/len;
		y=y/len;
		z=z/len;
	}
	vec cp(vec v)
	{
		return vec(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);
	}
	double dp(vec v)
	{
		return x*v.x+y*v.y+z*v.z;
	}
	vec add(vec v)
	{
		return vec(v.x+x,v.y+y,v.z+z);
	}
	vec sub(vec v)
	{
		return vec(-v.x+x,-v.y+y,-v.z+z);
	}
	vec scale(double s)
	{
		return vec(s*x,s*y,s*z);
	}

	vec rotateRodrigues(double angle,vec axis,vec v)
	{
	    angle*=(pi/180.0);
	    vec ret=v.scale(cos(angle));
	    vec cross=axis.cp(v);
	    cross=cross.scale(sin(angle));
	    ret=ret.add(cross);
	    vec tmp=axis.scale(axis.dp(v));
	    tmp=tmp.scale(1-cos(angle));
	    ret=ret.add(tmp);
	    ret.norm();
	    return ret;
	}
};

class mat4x4
{
public:
    double m[4][4];
    mat4x4()
    {
        int i,j;
        for(i=0;i<4;i++)for(j=0;j<4;j++)m[i][j]=(i==j)?1.0:0.0;
    }
    hPoint transformPoint(hPoint point)
    {
        int i,j;
        double sum;
        hPoint ret;
        for(i=0;i<4;i++)
        {
            sum=0;
            for(j=0;j<4;j++)sum+=m[i][j]*point.p[j];
            ret.p[i]=sum;
        }
        ret.normalize();
        return ret;
    }
    mat4x4 multiply(mat4x4 mat)
    {
        int i,j,k;
        double sum;
        mat4x4 ret;
        for(i=0;i<4;i++)for(j=0;j<4;j++)
        {
            sum=0;
            for(k=0;k<4;k++)sum+=m[i][k]*mat.m[k][j];
            ret.m[i][j]=sum;
        }
        return ret;
    }
    void generateTranslationMatrix(double tx, double ty, double tz)
    {
        int i,j;
        for(i=0;i<4;i++)for(j=0;j<4;j++)m[i][j]=(i==j)?1.0:0.0;
        m[0][3]=tx;
        m[1][3]=ty;
        m[2][3]=tz;
    }
    void generateScaleMatrix(double sx, double sy, double sz)
    {
        int i,j;
        for(i=0;i<4;i++)for(j=0;j<4;j++)m[i][j]=(i==j)?1.0:0.0;
        m[0][0]=sx;
        m[1][1]=sy;
        m[2][2]=sz;
    }
    void generateRotationMatrix(double angle, double ax, double ay, double az)
    {
        vec axis=vec(ax,ay,az);
        axis.norm();
        vec newI=axis.rotateRodrigues(angle,axis,vec(1,0,0));
        vec newJ=axis.rotateRodrigues(angle,axis,vec(0,1,0));
        vec newK=axis.rotateRodrigues(angle,axis,vec(0,0,1));
        m[0][0]=newI.x;
        m[1][0]=newI.y;
        m[2][0]=newI.z;
        m[0][1]=newJ.x;
        m[1][1]=newJ.y;
        m[2][1]=newJ.z;
        m[0][2]=newK.x;
        m[1][2]=newK.y;
        m[2][2]=newK.z;
    }
    void print()
    {
        int i,j;
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                printf("%.3lf\t",m[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
};

point eye,look;
vec up;
vec l,r,u;
double near,far;
double fovx,fovy,aspectratio;
stack<mat4x4> S;
stack<int> pushIndices;
vector<hPoint> stage1;
vector<hPoint> stage2;
vector<hPoint> stage3;
FILE* fin;
FILE* s1;
FILE* s2;
FILE* s3;

void model()
{
	int i;
    char command[20];
    S.push(mat4x4());
    while(true)
    {
        fscanf(fin," %s",command);
        if(strcmp("triangle",command)==0)
        {
            int l=3;
            while(l--)
            {
                double x,y,z;
                fscanf(fin," %lf %lf %lf",&x,&y,&z);
                hPoint p(x,y,z);
                mat4x4 mat=S.top();
                p=mat.transformPoint(p);
                stage1.push_back(p);
            }
        }
        else if(strcmp("translate",command)==0)
        {
            double tx,ty,tz;
            fscanf(fin," %lf %lf %lf",&tx,&ty,&tz);
            mat4x4 tmp=mat4x4();
            tmp.generateTranslationMatrix(tx,ty,tz);
            mat4x4 top=S.top();
            tmp=top.multiply(tmp);
            S.push(tmp);
        }
        else if(strcmp("scale",command)==0)
        {
            double sx,sy,sz;
            fscanf(fin," %lf %lf %lf",&sx,&sy,&sz);
            mat4x4 tmp=mat4x4();
            tmp.generateScaleMatrix(sx,sy,sz);
            mat4x4 top=S.top();
            tmp=top.multiply(tmp);
            S.push(tmp);
        }
        else if(strcmp("rotate",command)==0)
        {
            double angle,ax,ay,az;
            fscanf(fin," %lf %lf %lf %lf",&angle,&ax,&ay,&az);
            mat4x4 tmp=mat4x4();
            tmp.generateRotationMatrix(angle,ax,ay,az);
            //tmp.print();
            mat4x4 top=S.top();
            tmp=top.multiply(tmp);
            S.push(tmp);
        }
        else if(strcmp("push",command)==0)
        {
            int s=S.size();
            pushIndices.push(s);
        }
        else if(strcmp("pop",command)==0)
        {
            if(pushIndices.size()==0)
            {
                printf("error\n");
                continue;
            }
            int s=S.size();
            int restore=pushIndices.top();
            pushIndices.pop();
            int k=s-restore;
            while(k--)S.pop();
        }
        else if(strcmp("end",command)==0)break;
    }
	for(i=0;i<stage1.size();i++)
	{
		fprintf(s1,"%.7lf %.7lf %.7lf\n",stage1[i].p[0],stage1[i].p[1],stage1[i].p[2]);
		if(i%3==2)fprintf(s1,"\n");
	}
}

void view()
{
    l=vec(look.x-eye.x,look.y-eye.y,look.z-eye.z);
	l.norm();
	r=l.cp(up);
    r.norm();
    u=r.cp(l);
    u.norm();
    mat4x4 t=mat4x4();
    t.generateTranslationMatrix(-eye.x,-eye.y,-eye.z);
    mat4x4 view=mat4x4();
    view.m[2][0]=-l.x;
    view.m[2][1]=-l.y;
    view.m[2][2]=-l.z;
    view.m[0][0]=r.x;
    view.m[0][1]=r.y;
    view.m[0][2]=r.z;
    view.m[1][0]=u.x;
    view.m[1][1]=u.y;
    view.m[1][2]=u.z;
    view=view.multiply(t);
    int i;
    for(i=0;i<stage1.size();i++)
    {
        stage2.push_back(view.transformPoint(stage1[i]));
    }
	for(i=0;i<stage2.size();i++)
	{
		fprintf(s2,"%.7lf %.7lf %.7lf\n",stage2[i].p[0],stage2[i].p[1],stage2[i].p[2]);
		if(i%3==2)fprintf(s2,"\n");
	}
    return;
}

void projection()
{
	fovx=fovy*aspectratio;
    fovx*=(pi/180.0);
    fovy*=(pi/180.0);
    double t=near*tan(fovy/2.0);
    double r=near*tan(fovx/2.0);

    mat4x4 projection=mat4x4();
    projection.m[0][0]=near/r;
    projection.m[1][1]=near/t;
    projection.m[2][2]=(-1.0*(far+near))/(far-near);
    projection.m[2][3]=(-2.0*far*near)/(far-near);
    projection.m[3][3]=0.0;
    projection.m[3][2]=-1.0;

    int i;
    for(i=0;i<stage2.size();i++)
    {
        stage3.push_back(projection.transformPoint(stage2[i]));
    }
	for(i=0;i<stage3.size();i++)
	{
		fprintf(s3,"%.7lf %.7lf %.7lf\n",stage3[i].p[0],stage3[i].p[1],stage3[i].p[2]);
		if(i%3==2)fprintf(s3,"\n");
	}
    return;
}

void init()
{
	fin=fopen("scene.txt","r");

	s1=fopen("stage1.txt","w");
	s2=fopen("stage2.txt","w");
	s3=fopen("stage3.txt","w");

	if(fin==NULL||s1==NULL||s2==NULL||s3==NULL)
		printf("Error opening file");

	fscanf(fin," %lf %lf %lf",&eye.x,&eye.y,&eye.z);
	fscanf(fin," %lf %lf %lf",&look.x,&look.y,&look.z);
	fscanf(fin," %lf %lf %lf",&up.x,&up.y,&up.z);
	fscanf(fin," %lf %lf",&fovy,&aspectratio);
	fscanf(fin," %lf %lf",&near,&far);

	stage1.clear();
	stage2.clear();
	stage3.clear();
}

void finish()
{
    stage1.clear();
	stage2.clear();
	stage3.clear();

	while(!S.empty())S.pop();

	fclose(fin);
	fclose(s1);
	fclose(s2);
	fclose(s3);
}

int main()
{
	init();

	model();
	view();
	projection();

	finish();

	return 0;
}
