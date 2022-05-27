#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <GL/gl.h>
#include <GL/glut.h>
#include <vector>
#include <math.h>
#include "bitmap_image.hpp"
//extern int recursion_level;
class Vectors
{
public:
	double x,y,z;
	Vectors(double x,double y,double z){
		this->x = x;
		this->y = y;
		this->z = z;
	}
	Vectors(){
	}
	Vectors operator + (Vectors v) {
        Vectors V(x + v.x, y + v.y, z + v.z);
        return V;
    }

    Vectors operator - (Vectors v) {
        Vectors V(x - v.x, y - v.y, z - v.z);
        return V;
    }

    Vectors operator * (double d) {
        Vectors V(x*d, y*d, z*d);
        return V;
    }

    Vectors operator / (double d) {
        Vectors V(x/d, y/d, z/d);
        return V;
    }
};
void normalizer(Vectors * v)
{
    double val = sqrt((v->x*v->x)+(v->y*v->y)+(v->z*v->z));
    v->x = v->x /val;
    v->y = v->y /val;
    v->z = v->z /val;
}
Vectors crossProduct(Vectors v1,Vectors v2)
{
    Vectors v3(0,0,0);
    v3.x = v1.y * v2.z - v1.z * v2.y;
    v3.y = v1.z * v2.x - v1.x * v2.z;
    v3.z = v1.x * v2.y - v1.y * v2.x;
    return v3;
}
double dotProduct(Vectors v1,Vectors v2)
{
    double d;
    d = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
    return d;
}
class Rey{
public:
    Vectors start;
    Vectors direction;
    Rey(Vectors s,Vectors d){
        start = s;
        direction = d;
        //printf("Yay Rey created");
        normalizer(&direction);
    }
};
class Object{
public:
	Vectors reference_point;
	double height, width, length;
	double source_factor = 1.0, refraction_idx = 1.5;
	int Shine;
	double color[3];
	double co_efficients[4];
	Object(){
	}
	virtual void draw(){
	}
	virtual double intersectT(Rey r){
	    return -1;
	}
	virtual double intersect(Rey r, double current_color[3], int level){
	    return -1;
	}
	virtual Vectors normal(Vectors intersection){
	    Vectors v(0,0,0);
	    return v;
	}
	void setColor(double i, double j, double k){
        color[0] = i; color[1] = j; color[2] = k;
	}
	void setShine(int s){
        Shine = s;
	}
	void setCoEfficients(double a,double d,double s,double r){
        co_efficients[0] = a; co_efficients[1] = d;
        co_efficients[2] = s; co_efficients[3] = r;
	}
	Vectors reflection(Rey r,Vectors n){
        //std::cout<<dotProduct(r->direction,n)<<std::endl;
        double d = dotProduct(r.direction,n);
        Vectors refl;
        refl = r.direction - n * 2.0 *  d;
        normalizer(&refl);
        return refl;
	}
	Vectors refraction(Rey r,Vectors n){
        Vectors refr(0.0,0.0,0.0);
        double d = dotProduct(n,r.direction);
        double val = 1.0 - refraction_idx * refraction_idx * (1.0 - d * d);
        if(val >= 0){
        	refr = r.direction * refraction_idx - n * (refraction_idx * d + sqrt(val));
            normalizer(&refr);
        }
        return refr;
	}
};
std::vector <Vectors> lights;
std::vector <Object *> objects;
class Sphere: public Object{
public:
	Sphere(Vectors Center,double Radius){
		reference_point.x = Center.x;
		reference_point.y = Center.y;
		reference_point.z = Center.z;
		length=Radius;
	}
	void draw(){
	    //printf("entered here\n");
	    glColor3f(color[0], color[1], color[2]);
        Vectors points[100][100];
        double h, r;
        int slices = 32, stacks = 32;
        //generate points
        for(int i=0; i<=stacks; i++) {
            h = length * sin(((double)i/(double)stacks)*(pi/2));
            r = length * cos(((double)i/(double)stacks)*(pi/2));
            for(int j=0; j<=slices; j++) {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(int i=0; i<stacks; i++) {
            for(int j=0; j<slices; j++) {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x+reference_point.x,points[i][j].y+reference_point.y,points[i][j].z+reference_point.z);
                    glVertex3f(points[i][j+1].x+reference_point.x,points[i][j+1].y+reference_point.y,points[i][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j+1].x+reference_point.x,points[i+1][j+1].y+reference_point.y,points[i+1][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j].x+reference_point.x,points[i+1][j].y+reference_point.y,points[i+1][j].z+reference_point.z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x+reference_point.x,points[i][j].y+reference_point.y,-points[i][j].z+reference_point.z);
                    glVertex3f(points[i][j+1].x+reference_point.x,points[i][j+1].y+reference_point.y,-points[i][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j+1].x+reference_point.x,points[i+1][j+1].y+reference_point.y,-points[i+1][j+1].z+reference_point.z);
                    glVertex3f(points[i+1][j].x+reference_point.x,points[i+1][j].y+reference_point.y,-points[i+1][j].z+reference_point.z);
                }glEnd();
            }
        }
	}
    double intersectT(Rey r){
        Vectors starting_point;
        starting_point = r.start - this->reference_point;
        double A = dotProduct(r.direction,r.direction);
        double B = 2.0* dotProduct(r.direction,starting_point);
        double C = dotProduct(starting_point,starting_point) - this->length*this->length;
        double D = B*B - 4.0*A*C;
        if(D < 0) return -1;
        double root1 = (-B + sqrt(D))/(2.0*A);
        double root2 = (-B - sqrt(D))/(2.0*A);
        if(root1 < root2) return root1;
        else return root2;
    }
    Vectors normal(Vectors intersection){
        Vectors n;
        n = intersection - reference_point;
        //normalizer(&n);
        //std::cout<<n.x<<" "<<n.y<<" "<<n.z<<std::endl;
        return n;
    };
    double intersect(Rey r, double c[3], int level){
        //r is viewer's ray
        double t = intersectT(r);
        if(t <= 0){
            return -1;
        }
        if(level == 0){
            return t;
        }
        //std::cout<<level<<std::endl;
        for(int k=0;k<3;k++){
            c[k] = source_factor * color[k] * co_efficients[0];
        }
        Vectors intersecting_point;
        intersecting_point = r.start + r.direction * t;
        Vectors n = normal(intersecting_point);
        double value = dotProduct(r.direction,n);
        if (value > 0){
            n = n * (-1.0);
        }
        //Vectors refl = reflection(r,n);
        Vectors refl;Vectors refr;
        //normalizer(&refl);
        //std::cout<<"Normal "<<n.x<<" "<<n.y<<" "<<n.z<<std::endl;
        for(int i=0;i<lights.size();i++){
            Vectors dir = lights[i] - intersecting_point;
            double len = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
            normalizer(&dir);
            Vectors start = intersecting_point + dir * 1.0;
            Rey L(start,dir);
            //L is light source's ray
            bool obscure = true;
            for(int j=0;j<objects.size();j++){//check if obscured
                double t_prime = objects[j]->intersectT(L);
                if(t_prime > 0 && fabs(t_prime) < len) {
                    obscure = false;
                    break;
                }
                //obscure = true;
                //break;
            }
            refl = reflection(r,n);
            //normalizer(&refl);
            //refl = reflection(L,n);
            //refr = refraction(L,n);
            //refr = refraction(r,n);
            //normalizer(&refr);
            if(obscure){
                double lambart = std::max(dotProduct(L.direction,n),0.0);
                double phong = std::max(pow(dotProduct(refl,r.direction),Shine),0.0);
                //std::cout<<"Lambart and Phong "<<lambart<<" "<<phong<<std::endl;
                for(int k=0;k<3;k++){
                    c[k] += (lambart * co_efficients[1]* color[k]);
                    c[k] += (phong * co_efficients[2] * color[k]);
                }
            }
        }
        if(level < recursion_level){
            //Reflection
            Vectors start = intersecting_point + refl * 1.0;
            Rey reflected_rey(start,refl);
            int nearestL=-1;
            double minvalL = 9999;
            double reflected_color[3];
            for(int k=0;k<objects.size();k++){
                double t_ref = objects[k]->intersectT(reflected_rey);
                if(t_ref <= 0){
                    continue;
                }
                else if(t_ref < minvalL){
                    minvalL = t_ref;
                    nearestL = k;
                }
            }
            if(nearestL != -1){
                objects[nearestL]->intersect(reflected_rey,reflected_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (reflected_color[k] * co_efficients[3]);
                }
            }
            //Refraction
            refr = refraction(r,n);
            start = intersecting_point + refr * 1.0;
            Rey refracted_rey(start,refr);
            int nearestR=-1;
            double minvalR = 9999;
            double refracted_color[3];
            for(int k=0;k<objects.size();k++){
                double t_raf = objects[k]->intersectT(refracted_rey);
                if(t_raf <= 0){
                    continue;
                }
                else if(t_raf < minvalR){
                    minvalR = t_raf;
                    nearestR = k;
                }
            }
            if(nearestR != -1){
                objects[nearestR]->intersect(refracted_rey,refracted_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (refracted_color[k] * refraction_idx);
                }
            }
        }
        for(int k=0;k<3;k++){
            if (c[k] > 1){
                c[k] = 1;
            }
            if(c[k] < 0){
                c[k] = 0;
            }
        }
        return t;
	}
};
class Triangle: public Object{
public:
    Vectors A;Vectors B;Vectors C;
    Triangle(Vectors a,Vectors b,Vectors c){
        A = a;
        B = b;
        C = c;
    }
    void draw(){
        //printf("entered here\n");
	    glColor3f(color[0],color[1],color[2]);
	    glPushMatrix();{
            glBegin(GL_TRIANGLES);{
                glVertex3f(A.x, A.y, A.z);
                glVertex3f(B.x, B.y, B.z);
                glVertex3f(C.x, C.y, C.z);
            }glEnd();
	    }glPopMatrix();
    }
    Vectors normal(Vectors intersection){
        Vectors p = B-A;
        Vectors q = C-A;
        Vectors n = crossProduct(p,q);
        normalizer(&n);
        return n;
    }
    double intersectT(Rey r){
        const double EPSILON = 0.0000001;
        double a,f,u,v;
        Vectors edge1 = B - A;
        Vectors edge2 = C - A;
        Vectors edge = crossProduct(r.direction,edge2);
        a = dotProduct(edge1,edge);
        if (a > -EPSILON && a < EPSILON)
            return -1;
        f = 1/a;
        Vectors e = r.start - A;
        u = f * dotProduct(e,edge);
        if (u < 0.0 || u > 1.0)
            return -1;
        Vectors q = crossProduct(e,edge1);
        v = f * dotProduct(r.direction,q);
        if (v < 0.0 || u + v > 1.0)
            return -1;
        // At this stage we can compute t to find out where the intersection point is on the line.
        double t = f * dotProduct(edge2,q);
        if (t > EPSILON) // ray intersection
        {
            return t;
        }
        else{ // This means that there is a line intersection but not a ray intersection.
            return -1;
        }
    }
    double intersect(Rey r, double c[3], int level){
        //r is viewer's ray
        double t = intersectT(r);
        if(t <= 0){
            return -1;
        }
        if(level == 0){
            return t;
        }
        //std::cout<<level<<std::endl;
        for(int k=0;k<3;k++){
            c[k] = color[k] * co_efficients[0];
        }
        Vectors intersecting_point;
        intersecting_point = r.start + r.direction * t;
        Vectors n = normal(intersecting_point);
        double value = dotProduct(r.direction,n);
        if (value > 0){
            n = n * (-1.0);
        }
        Vectors refl;Vectors refr;
        //Vectors refl = reflection(r,n);
        //normalizer(&refl);
        //std::cout<<"Normal "<<n.x<<" "<<n.y<<" "<<n.z<<std::endl;
        for(int i=0;i<lights.size();i++){
            Vectors dir = lights[i] - intersecting_point;
            double len = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
            normalizer(&dir);
            /**double value = dotProduct(dir,n);
            if(value>0){
                n = n * (-1.0);
            }**/
            Vectors start = intersecting_point + dir * 1.0;
            Rey L(start,dir);
            //L is light source's ray
            bool obscure = true;
            for(int j=0;j<objects.size();j++){//check if obscured
                double t_prime = objects[j]->intersectT(L);
                if(t_prime > 0 && fabs(t_prime) < len) {
                    obscure = false;
                    break;
                }
                //obscure = true;
                //break;
            }
            //refl = reflection(L,n);
            refl = reflection(r,n);
            //refr = refraction(L,n);
            //refr = refraction(r,n);
            if(obscure){
                double lambart = std::max(dotProduct(L.direction,n),0.0);
                double phong = std::max(pow(dotProduct(refl,r.direction),Shine),0.0);
                //std::cout<<"Lambart and Phong "<<lambart<<" "<<phong<<std::endl;
                for(int k=0;k<3;k++){
                    c[k] += (lambart * co_efficients[1]* color[k]);
                    c[k] += (phong * co_efficients[2] * color[k]);
                }
            }
        }
        if(level < recursion_level){
            //Reflection
            Vectors start = intersecting_point + refl * 1.0;
            Rey reflected_rey(start,refl);
            int nearestL=-1;
            double minvalL = 9999;
            double reflected_color[3];
            for(int k=0;k<objects.size();k++){
                double t_ref = objects[k]->intersectT(reflected_rey);
                if(t_ref <= 0){
                    continue;
                }
                else if(t_ref < minvalL){
                    minvalL = t_ref;
                    nearestL = k;
                }
            }
            if(nearestL != -1){
                objects[nearestL]->intersect(reflected_rey,reflected_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (reflected_color[k] * co_efficients[3]);
                }
            }
            //Refraction
            /**refr = refraction(r,n);
            start = intersecting_point + refr * 1.0;
            Rey refracted_rey(start,refr);
            int nearestR=-1;
            double minvalR = 9999;
            double refracted_color[3];
            for(int k=0;k<objects.size();k++){
                double t_raf = objects[k]->intersectT(refracted_rey);
                if(t_raf <= 0){
                    continue;
                }
                else if(t_raf < minvalR){
                    minvalR = t_raf;
                    nearestR = k;
                }
            }
            if(nearestR != -1){
                objects[nearestR]->intersect(refracted_rey,refracted_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (refracted_color[k] * refraction_idx);
                }
            }**/
        }
        for(int k=0;k<3;k++){
            if (c[k] > 1){
                c[k] = 1;
            }
            if(c[k] < 0){
                c[k] = 0;
            }
        }
        return t;
    }
};
class Floor: public Object{
public:
    bitmap_image shades;
    double img_height,img_width;
    Floor(double FloorWidth,double TileWidth){
        reference_point.x = -FloorWidth/2;
        reference_point.y = -FloorWidth/2;
        reference_point.z = 0;
        length=TileWidth;
        shades = bitmap_image("D:\\Study\\GraphicsCodes\\RayTracer\\shades.bmp");
        img_height = shades.height()/1000.0;
        img_width = shades.width()/1000.0;
    }
    void draw(){
        int GridSize = fabs(reference_point.x *2 /length);
        for(int x =0;x<GridSize;x++){
            for(int y =0;y<GridSize;y++){
                if ((x+y)%2 == 0){ //modulo 2
                    //printf("yes\n");
                    glColor3f(1.0f,1.0f,1.0f); //white
                }
                else{
                    glColor3f(0.0f,0.0f,0.0f); //black
                    //printf("no\n");
                }
                glBegin(GL_QUADS);
                {
                    glVertex3f(reference_point.x+length*x,reference_point.y+length*y, reference_point.z);
                    glVertex3f(reference_point.x+length*(x+1),reference_point.y+length*y, reference_point.z);
                    glVertex3f(reference_point.x+length*(x+1), reference_point.y+length*(y+1), reference_point.z);
                    glVertex3f(reference_point.x+length*x, reference_point.y+length*(y+1), reference_point.z);
                }
                glEnd();
            }
        }
    }
    Vectors normal(Vectors intersection){
        Vectors n(0,0,1);
        normalizer(&n);
        return n;
    }
    double intersectT(Rey r){
        Vectors n = normal(reference_point);
        double t = dotProduct(n,r.start) * (-1.0) / dotProduct(n,r.direction);
        return t;
    }
    double intersect(Rey r, double c[3], int level){
        double t = intersectT(r);
        if(t <= 0) return -1;
        if(level == 0) return t;
        Vectors intersecting_point = r.start + r.direction * t;
        /**double x_min = reference_point.x;
        double x_max = x_min * (-1.0);
        double y_min = reference_point.y;
        double y_max = y_min * (-1.0);
        if(x_min > intersecting_point.x || x_max < intersecting_point.x ||
                y_min > intersecting_point.y || y_max < intersecting_point.y) {
            return -1;
        }**/
        int x = (intersecting_point.x-reference_point.x) /this->length;
        int y = (intersecting_point.y-reference_point.y) / this->length;
        if((x+y)%2 == 0){
            color[0] = color[1] = color[2] = 0;
        }
        else{
            color[0] = color[1] = color[2] = 1;
        }
        unsigned char red,green,blue;
        shades.get_pixel(x, y, red, green, blue);
        double rgb[] = {red, green, blue};
        for(int k=0;k<3;k++){
            c[k] = color[k] * co_efficients[0] * (rgb[k] / 255.0);
        }
        //Illumination
        Vectors n = normal(intersecting_point);
        double value = dotProduct(r.direction,n);
        if (value > 0){
            n = n * (-1.0);
        }
        //Vectors refl = reflection(r,n);
        //normalizer(&refl);
        Vectors refl;Vectors refr;
        //std::cout<<"Normal "<<n.x<<" "<<n.y<<" "<<n.z<<std::endl;
        for(int i=0;i<lights.size();i++){
            Vectors dir = lights[i] - intersecting_point;
            double len = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
            normalizer(&dir);
            Vectors start = intersecting_point + dir * 1.0;
            Rey L(start,dir);
            //L is light source's ray
            bool obscure = true;
            for(int j=0;j<objects.size();j++){//check if obscured
                double t_prime = objects[j]->intersectT(L);
                if(t_prime > 0 && fabs(t_prime) <= len) {
                    obscure = false;
                    break;
                }
                //obscure = true;
                //break;
            }
            //refl = reflection(L,n);
            //refr = refraction(L,n);
            refl = reflection(r,n);
            //refr = refraction(r,n);
            if(obscure){
                double lambart = std::max(dotProduct(L.direction,n),0.0);
                double phong = std::max(pow(dotProduct(refl,r.direction),Shine),0.0);
                //std::cout<<"Lambart and Phong "<<lambart<<" "<<phong<<std::endl;
                for(int k=0;k<3;k++){
                    c[k] += (lambart * co_efficients[1]* color[k]);
                    c[k] += (phong * co_efficients[2] * color[k]);
                }
            }
        }
        if(level < recursion_level){
            //Reflection
            Vectors start = intersecting_point + refl * 1.0;
            Rey reflected_rey(start,refl);
            int nearestL=-1;
            double minvalL = 9999;
            double reflected_color[3];
            for(int k=0;k<objects.size();k++){
                double t_ref = objects[k]->intersectT(reflected_rey);
                if(t_ref <= 0){
                    continue;
                }
                else if(t_ref < minvalL){
                    minvalL = t_ref;
                    nearestL = k;
                }
            }
            if(nearestL != -1){
                objects[nearestL]->intersect(reflected_rey,reflected_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (reflected_color[k] * co_efficients[3]);
                }
            }
            //Refraction
            /**refr = refraction(r,n);
            start = intersecting_point + refr * 1.0;
            Rey refracted_rey(start,refr);
            int nearestR=-1;
            double minvalR = 9999;
            double refracted_color[3];
            for(int k=0;k<objects.size();k++){
                double t_raf = objects[k]->intersectT(refracted_rey);
                if(t_raf <= 0){
                    continue;
                }
                else if(t_raf < minvalR){
                    minvalR = t_raf;
                    nearestR = k;
                }
            }
            if(nearestR != -1){
                objects[nearestR]->intersect(refracted_rey,refracted_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (refracted_color[k] * refraction_idx);
                }
            }**/
        }
        for(int k=0;k<3;k++){
            if (c[k] > 1){
                c[k] = 1;
            }
            if(c[k] < 0){
                c[k] = 0;
            }
        }
        return t;
    }
};
class GeneralQuad:public Object{
public:
    double A,B,C,D,E,F,G,H,I,J;
    GeneralQuad(double s[10],Vectors v,double a,double b,double c){
        A = s[0];
        B = s[1];
        C = s[2];
        D = s[3];
        E = s[4];
        F = s[5];
        G = s[6];
        H = s[7];
        I = s[8];
        J = s[9];
        reference_point = v;
        length = a;
        width = b;
        height = c;
    }
    void draw(){}
    Vectors normal(Vectors v){
        double a = 2.0 *A*v.x + D*v.y + F*v.z  + G;
        double b = 2.0*B*v.y + D*v.x + E*v.z  + H;
        double c = 2.0*C*v.z + E*v.y + F*v.x  + I;
        Vectors n(a,b,c);
        normalizer(&n);
    }
    double intersectT(Rey r){
        double a,b,c;
        //squared terms
        a = A*r.direction.x*r.direction.x + B*r.direction.y*r.direction.y + C*r.direction.z*r.direction.z;
        b = 2.0*(A*r.start.x*r.direction.x + B*r.start.y*r.direction.y + C*r.start.z*r.direction.z);
        c = A*r.start.x*r.start.x + B*r.start.y*r.start.y + C*r.start.z*r.start.z;
        //multiply terms
        a += D*r.direction.x*r.direction.y + E*r.direction.y*r.direction.z + F*r.direction.z*r.direction.x;
        b += D*(r.start.x*r.direction.y + r.direction.x*r.start.y)
            + E*(r.start.y*r.direction.z + r.direction.y*r.start.z)
            + F*(r.start.z*r.direction.x + r.direction.z*r.start.x);
        c += D*r.start.x*r.start.y + E*r.start.y*r.start.z + F*r.start.z*r.start.x;
        //constant terms
        b += G*r.direction.x + H*r.direction.y + I*r.direction.z;
        c += G*r.start.x + H*r.start.y + I*r.start.z + J;
        double d = b*b - 4.0*a*c;
        if(d < 0){
            return -1;
        }
        double t1 = (-b + sqrt(d))/(2.0*a);
        double t2 = (-b - sqrt(d))/(2.0*a);
        Vectors Point1 = r.start + r.direction*t1;
        Vectors Point2 = r.start + r.direction*t2;
        double x_min = reference_point.x;
        double x_max = x_min + length;
        double y_min = reference_point.y;
        double y_max = y_min + width;
        double z_min = reference_point.z;
        double z_max = z_min + height;
        bool flag1 = (length > 0 && (x_min > Point1.x || x_max < Point1.x) || width > 0 && (y_min > Point1.y || y_max < Point1.y) || height > 0 && (z_min > Point1.z || z_max < Point1.z));
        bool flag2 = (length > 0 && (x_min > Point2.x || x_max < Point2.x) || width > 0 && (y_min > Point2.y || y_max < Point2.y) || height > 0 && (z_min > Point2.z || z_max < Point2.z));
        if(flag1 && flag2){
            return -1;
        }
        else if(flag1){
            return t2;
        }
        else if(flag2){
            return t1;
        }
        else{
            if(t2 > t1) return t1;
            else return t2;
        }
    }
    double intersect(Rey r, double c[3], int level){
        //r is viewer's ray
        double t = intersectT(r);
        if(t <= 0){
            return -1;
        }
        if(level == 0){
            return t;
        }
        //std::cout<<level<<std::endl;
        for(int k=0;k<3;k++){
            c[k] = color[k] * co_efficients[0];
        }
        Vectors intersecting_point;
        intersecting_point = r.start + r.direction * t;
        Vectors n = normal(intersecting_point);
        //normalizer(&n);
        double value = dotProduct(r.direction,n);
        if (value > 0){
            n = n * (-1.0);
        }
        //Vectors refl = reflection(r,n);
        Vectors refl;Vectors refr;
        /**double value = dotProduct(r.direction,n);
        if(value > 0){
            n = n * (-1.0);
        }**/
        //normalizer(&refl);
        //std::cout<<"Normal "<<n.x<<" "<<n.y<<" "<<n.z<<std::endl;
        for(int i=0;i<lights.size();i++){
            Vectors dir = lights[i] - intersecting_point;
            double len = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
            normalizer(&dir);
            Vectors start = intersecting_point + dir * 1.0;
            Rey L(start,dir);
            //L is light source's ray
            bool obscure = true;
            for(int j=0;j<objects.size();j++){//check if obscured
                double t_prime = objects[j]->intersectT(L);
                if(t_prime > 0 && fabs(t_prime) <= len) {
                    obscure = false;
                    break;
                }
                //obscure = true;
                //break;
            }
            refl = reflection(r,n);
            //normalizer(&refl);
            //refl = reflection(L,n);
            //refr = refraction(L,n);
            //refr = refraction(r,n);
            //normalizer(&refr);
            if(obscure){
                double lambart = std::max(dotProduct(L.direction,n),0.0);
                double phong = std::max(pow(dotProduct(refl,r.direction),Shine),0.0);
                //std::cout<<"Lambart and Phong "<<lambart<<" "<<phong<<std::endl;
                for(int k=0;k<3;k++){
                    c[k] += (lambart * co_efficients[1]* color[k]);
                    c[k] += (phong * co_efficients[2] * color[k]);
                }
            }
        }
        if(level < recursion_level){
            //Reflection
            Vectors start = intersecting_point + refl * 1.0;
            Rey reflected_rey(start,refl);
            int nearestL=-1;
            double minvalL = 9999;
            double reflected_color[3];
            for(int k=0;k<objects.size();k++){
                double t_ref = objects[k]->intersectT(reflected_rey);
                if(t_ref <= 0){
                    continue;
                }
                else if(t_ref < minvalL){
                    minvalL = t_ref;
                    nearestL = k;
                }
            }
            if(nearestL != -1){
                objects[nearestL]->intersect(reflected_rey,reflected_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (reflected_color[k] * co_efficients[3]);
                }
            }
            //Refraction
            /**refr = refraction(r,n);
            start = intersecting_point + refr * 1.0;
            Rey refracted_rey(start,refr);
            int nearestR=-1;
            double minvalR = 9999;
            double refracted_color[3];
            for(int k=0;k<objects.size();k++){
                double t_raf = objects[k]->intersectT(refracted_rey);
                if(t_raf <= 0){
                    continue;
                }
                else if(t_raf < minvalR){
                    minvalR = t_raf;
                    nearestR = k;
                }
            }
            if(nearestR != -1){
                objects[nearestR]->intersect(refracted_rey,refracted_color,level+1);
                for(int k=0;k<3;k++){
                    c[k] += (refracted_color[k] * refraction_idx);
                }
            }**/
        }
        for(int k=0;k<3;k++){
            if (c[k] > 1){
                c[k] = 1;
            }
            if(c[k] < 0){
                c[k] = 0;
            }
        }
        return t;
	}
};
