#include <windows.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#define pi (2*acos(0.0))
#define window_width 500
#define window_height 500
#define fovy 80
using namespace std;
double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double a,r;
double A;
int image_width;
int image_height;
int recursion_level;
int *** frame_buffer;
//extern vector <struct Vectors> lights;
//extern vector <Object *> objects;
#include "class.hpp"

Vectors pos;
Vectors sop;
Vectors look;
Vectors rights;
Vectors up;

void drawPoints(Vectors v)
{
    glColor3f(0.8,0.5,0);
    glBegin(GL_QUADS);
    {
        glVertex3f(v.x+5, v.y, v.z+5);
        glVertex3f(v.x+5, v.y, v.z-5);
        glVertex3f(v.x-5, v.y, v.z-5);
        glVertex3f(v.x-5, v.y, v.z+5);
    }
    glEnd();
}
void drawAxes(){
	if(true)
	{

		glBegin(GL_LINES);{
		    glColor3f(0, 0, 1);
			glVertex3f(0,0,0);//x axis blue
			glVertex3f(150,0,0);
			glColor3f(0, 1, 0);
			glVertex3f(0,0,0);//y axis red
			glVertex3f(0,150,0);
		    glColor3f(1,0,0);
			glVertex3f(0,0,0);//z axis
			glVertex3f(0,0,150);
		}glEnd();
	}
}

void drawGrid(){
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 1, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes
                //XY Plane
				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
	if(drawgrid==2)
	{
		glColor3f(1, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes
                //YZ plane
				//lines parallel to -axis
				glVertex3f(0, -90, i*10);
				glVertex3f(0,  90, i*10);

				//lines parallel to -axis
				glVertex3f(0, i*10, -90);
				glVertex3f(0, i*10, 90);
			}
		}glEnd();
	}
	if(drawgrid==3)
	{
		glColor3f(0.6, 0.6, 1);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes
                //ZX Plane
				//lines parallel to -axis
				glVertex3f(-90, 0, i*10);
				glVertex3f(90,  0, i*10);

				//lines parallel to -axis
				glVertex3f(i*10, 0, -90);
				glVertex3f(i*10, 0, 90);
			}
		}glEnd();
	}
}

void draw_image(string name)
{
    bitmap_image image(image_width,image_height);//(width,height)
    //printf("drawing image\n");
    for(int i=0;i<image_width;i++){
        for(int j=0;j<image_height;j++){
            image.set_pixel(j,i,frame_buffer[i][j][0],frame_buffer[i][j][1],frame_buffer[i][j][2]);
        }
    }
    //cout<<name<<endl;
    image.save_image(name);
}

void capture(){
    //printf("hidden surface removal\n");
    frame_buffer = new int**[image_width];
    for(int i = 0; i < image_width; i++){
        frame_buffer[i] = new int*[image_height];
        for(int j = 0; j < image_height; j++){
            frame_buffer[i][j] = new int[3];
        }
    }
    for(int i = 0; i < image_width; i++){
        for(int j = 0; j < image_height; j++){
            for(int k = 0; k < 3; k++){
                frame_buffer[i][j][k] = 0;
            }
        }
    }
    double plane_distance = (window_height/2.0)/tan(fovy * pi/360);
    //printf("%lf",plane_distance);
    //cout<<"Plane Distance "<<plane_distance<<endl;
    Vectors toplook;
    toplook.x = pos.x + (look.x * plane_distance - rights.x * (window_width/2.0) + up.x * (window_height/2.0));
    toplook.y = pos.y + (look.y * plane_distance - rights.y * (window_width/2.0) + up.y * (window_height/2.0));
    toplook.z = pos.z + (look.z * plane_distance - rights.z * (window_width/2.0) + up.z * (window_height/2.0));
    //cout<<"Toplook "<<toplook.x<<" "<<toplook.y<<" "<<toplook.z<<endl;
    double du = (window_width*1.0)/image_width;
    double dv = (window_height*1.0)/image_height;
    for(int i=0; i<image_width;i++){
        for(int j=0; j<image_height;j++){
            Vectors corner;
            corner = toplook + rights * j * du - up * i * dv;
            Vectors dir; //corner - pos
            dir = corner - pos;
            Rey rey(pos,dir);
            int nearest=-1;
            double minval = 9999;
            double dummy_color[3];
            for (int k=0; k < objects.size(); k++){
                double t = objects[k]->intersect(rey,dummy_color,0);
                if(i%200==0 && j%200==0){
                    //cout<<"t is "<<t<<endl;
                }
                if(t <= 0){
                    continue;
                }
                else if(t < minval){
                    minval = t;
                    nearest = k;
                }
            }
            if(nearest != -1){
                //printf("no its not -1\n");
                double t = objects[nearest]->intersect(rey,dummy_color,1);
                frame_buffer[i][j][0] = dummy_color[0] * 255;
                frame_buffer[i][j][1] = dummy_color[1] * 255;
                frame_buffer[i][j][2] = dummy_color[2] * 255;
            }
        }
    }
    draw_image("output.bmp");
    printf("Generated Image\n");
}

void keyboardListener(unsigned char key, int x,int y){
    Vectors r1,l1,r2,l2,u1,u2,m,n;
	switch(key){
        case '0':
            capture();
            break;
	    case '1':
	        A = 0.05;
            r1 = rights;
            l1 = look;
            m = r1 * cos(-A);
            n = l1 * sin(-A);
            r2 = m + n;
            m = l1 * cos(-A);
            n = r1 * sin(-A);
            l2 = m - n;
            look = l2;
            rights = r2;
            break;
        case '2':
            A = 0.05;
            r1 = rights;
            l1 = look;
            m = r1* cos(A);
            n = l1* sin(A);
            r2 = m + n;
            m = l1* cos(A);
            n = r1* sin(A);
            l2 = m - n;
            look = l2;
            rights = r2;
            break;
        case '3':
            A = 0.05;
            u1 = up;
            l1 = look;
            m = l1* cos(A);
            n = u1 * sin(A);
            l2 = m + n;
            m = u1* cos(A);
            n = l1* sin(A);
            u2 = m - n;
            look = l2;
            up = u2;
            break;
        case '4':
            A = 0.05;
            u1 = up;
            l1 = look;
            m = l1* cos(-A);
            n = u1* sin(-A);
            l2 = m + n;
            m = u1* cos(-A);
            n = l1* sin(-A);
            u2 = m - n;
            look = l2;
            up = u2;
            break;
        case '5':
            A = 0.05;
            u1 = up;
            r1 = rights;
            m = u1* cos(A);
            n = r1* sin(A);
            u2 = m + n;
            m = r1* cos(A);
            n = u1* sin(A);
            r2 = m - n;
            up = u2;
            rights = r2;
            break;
        case '6':
            A = 0.05;
            u1 = up;
            r1 = rights;
            m = u1* cos(-A);
            n = r1* sin(-A);
            u2 = m + n;
            m = r1* cos(-A);
            n = u1* sin(-A);
            r2 = m - n;
            up = u2;
            rights = r2;
            break;
        case 'q':
            exit(0);
            break;
        case 'm':
            drawgrid=1;
            break;
        case 'n':
            drawgrid=2;
            break;
        case 'o':
            drawgrid=3;
        default:
            break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_UP:		// up arrow key
		    pos.x += look.x;
            pos.y += look.y;
            pos.z += look.z;
			break;
        case GLUT_KEY_DOWN:		//down arrow key
            pos.x -= look.x;
            pos.y -= look.y;
            pos.z -= look.z;
			break;
        case GLUT_KEY_LEFT:
            pos.x -= rights.x;
            pos.y -= rights.y;
            pos.z -= rights.z;
			break;
		case GLUT_KEY_RIGHT:
            pos.x += rights.x;
            pos.y += rights.y;
            pos.z += rights.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += up.x;
            pos.y += up.y;
            pos.z += up.z;
			break;

		case GLUT_KEY_PAGE_DOWN:
			pos.x -= up.x;
            pos.y -= up.y;
            pos.z -= up.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
		    if (a > 0){
                a = a - 0.5;
                r = r + 0.5;
		    }
			break;
		case GLUT_KEY_END:
		    if (r > 0){
                a = a + 0.5;
                r = r - 0.5;
		    }
			break;

		default:
			break;
	}
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?
    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(pos.x,pos.y,pos.z,0, 0, 0,up.x,up.y,up.z);
    gluLookAt(pos.x,pos.y,pos.z,pos.x + look.x, pos.y + look.y, pos.z + look.z,up.x,up.y,up.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	//drawGrid();
    for(int i=0; i< objects.size(); i++){
        Object * o;
        o = objects[i];
        o->draw();
    }
    for(int i=0; i< lights.size(); i++){
        Vectors v;
        v = lights[i];
        drawPoints(v);
    }
	glutSwapBuffers();
}

void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=0;
	cameraHeight=150.0;
	cameraAngle=1.0;
	a = 25;
	r = 5;
	angle=0;
	//image_height = 768;
	//image_width = 768;
    pos.x = 100;pos.y = 100;pos.z = 10;
    //pos.x =0; pos.y = -200; pos.z = 10;
    look.x = -1/sqrt(2);look.y = -1/sqrt(2);look.z = 0;
    rights.x = -1/sqrt(2);rights.y = 1/sqrt(2);rights.z = 0;
    up.x = 0;up.y = 0;up.z = 1;
	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(fovy,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void load_test_data(){
    //printf("Do Nothing");
    //cout<<"here"<<endl;
    Object *temp;
    struct Vectors Center;
    Center.x = 40; Center.y = 0; Center.z = 10;
    double Radius = 10;
    temp=new Sphere(Center, Radius); // Center(0,0,10), Radius 10
    temp->setColor(1,0,0);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);
    Center.x = -30; Center.y = 60; Center.z = 20;
    Radius = 20;
    temp=new Sphere(Center, Radius); // Center(0,0,10), Radius 10
    temp->setColor(0,1,0);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);
    Center.x = -15; Center.y = 15; Center.z = 15;
    Radius = 15;
    temp=new Sphere(Center, Radius); // Center(0,0,10), Radius 10
    temp->setColor(0,0,1);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);
    Vectors v1(50,30,0);Vectors v2(70,60,0);Vectors v3(50,45,50);
    temp = new Triangle(v1,v2,v3);
    temp->setColor(1,0,1);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);
    Vectors light1;
    light1.x = 70;light1.y = 70; light1.z = 70;
    lights.push_back(light1);
    Vectors light2;
    light2.x = -70;light2.y = 70; light2.z = 70;
    lights.push_back(light2);
    temp=new Floor(256, 8);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);
}

void load_actual_data(){
    string path;
    path = "D:\\Study\\GraphicsCodes\\RayTracer\\scene.txt";
    //ifstream infile(path.c_str());
    //string line;
    //std::getline(infile, line);
    //std::istringstream iss(line);
    freopen(path.c_str(),"r",stdin);
    cin >> recursion_level;
    cout << recursion_level << endl;
    //std::getline(infile, line);
    //std::istringstream iss(line);
    cin >> image_width;
    image_height = image_width;
    cout<<image_width<<endl;
    int object_count;
    //std::getline(infile, line);
    //std::istringstream iss(line);
    cin >> object_count;
        //cout<<object_count;
    string command;
    double a,b,c,Radius;
    Object * temp;
    for(int i=0;i<object_count;i++){
        //std::getline(infile, line);
        //std::istringstream iss(line);
        cin >> command;
        cout << command<<endl;
        if(command == "sphere"){
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >>c;
            Vectors Center(a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> Radius;
            temp = new Sphere(Center,Radius);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c;
            temp->setColor(a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c >> Radius;
            temp->setCoEfficients(a,b,c,Radius);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a;
            temp->setShine(a);
            objects.push_back(temp);
        }
        else if(command == "triangle"){
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c;
            Vectors A(a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c;
            Vectors B(a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a>> b >> c;
            Vectors C(a,b,c);
            temp = new Triangle(A,B,C);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c;
            temp->setColor(a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c >> Radius;
            temp->setCoEfficients(a,b,c,Radius);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a;
            temp->setShine(a);
            objects.push_back(temp);
        }
        else if(command == "general"){
            double s[10];
            for(int j=0; j<10; j++){
                //std::getline(infile, line);
                //std::istringstream iss(line);
                cin >> s[j];
            }
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c;
            Vectors v(a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c;
            temp = new GeneralQuad(s,v,a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c;
            temp->setColor(a,b,c);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a >> b >> c >> Radius;
            temp->setCoEfficients(a,b,c,Radius);
            //std::getline(infile, line);
            //std::istringstream iss(line);
            cin >> a;
            temp->setShine(a);
            objects.push_back(temp);
        }
    }
    //std::getline(infile, line);
    //std::istringstream iss(line);
    cin >> object_count;
    cout<<object_count<<endl;
    for(int i=0; i<object_count; i++){
        //std::getline(infile, line);
        //std::istringstream iss(line);
        cin >>a>>b>>c;
        Vectors light(a,b,c);
        lights.push_back(light);
    }
    fclose(stdin);
    //infile.close();
    //Object * temp;
    temp = new Floor(1000, 20);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);
}

void free_memory(){
    vector<Vectors>().swap(lights);
    vector<Object*>().swap(objects);
    //delete frame_buffer;
}

int main(int argc, char **argv){
    //load_test_data();
    load_actual_data();
	glutInit(&argc,argv);
	glutInitWindowSize(window_width,window_height);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
    glutCreateWindow("My OpenGL Program");
    init();
	glEnable(GL_DEPTH_TEST);	//enable Depth Testing
	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL
    free_memory();
	return 0;
}

