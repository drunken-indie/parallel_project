#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>
#include <ctime>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define MASTER 0    //master process
#define epsilon 0.000000000000000222
#define G 6.67408e-11 // The gravitational constant

int my_rank;

void Initialize_Input(int argc, char* argv[], int* numParticlesLight, int* numParticleMedium, int* numParticleHeavy, int* numSteps, int* subSteps, double* timeSubStep, int* width, int* height, char** output);
void Initialize_particles(body** bodies, int total_particles, int numParticlesLight, int numParticleMedium, int numParticleHeavy, int width, int height);
void print_status(body** bodies, int i);
void saveIt(body** bodies, char *filename, int width, int height, int total_particles);
void calculateForce(body** body1, body** body2, int a);
void updatePosVel(body** bodies, body** body1, int total_particles, double timeSubStep, int width, int height);
int main(int argc, char* argv[]){

    //variables
    int currentstep, currentsubstep;
    int width, height;
    int numParticlesLight;// = atoi(argv[1]);
    int numParticleMedium;// = atoi(argv[2]);
    int numParticleHeavy;// = atoi(argv[3]);

    int numSteps;// = atoi(argv[4]);
    int subSteps;// = atoi(argv[5]);
    double timeSubStep;//= atoi(argv[6]);
    char* output;
    int total_particles;
    body *bodies;
    char* filename;
    int filecounter=0;
    double start, end;
    double total_time=0;
    double fastest = 10000;
    double latest = 0;
    int i;
    

    
    //width=atoi(argv[7]);
    //height=atoi(argv[8]);

    
    
    //printf("sdsd \n"); 
    
    
    
    MPI_Init(&argc,&argv);

    int p;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    
    //create a tyoe for struct vec3
    const int nitems=3;
    int          blocklengths[3] = {1,1,1};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_vec3_type;
    MPI_Aint     offsets[3];

    offsets[0] = offsetof(vec3, x);
    offsets[1] = offsetof(vec3, y);
    offsets[2] = offsetof(vec3, z);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_vec3_type);
    MPI_Type_commit(&mpi_vec3_type);


    //create a tyoe for struct vec3
    const int nitems1=6;
    int          blocklengths1[6] = {1,1,1,1,1,1};
    MPI_Datatype types1[6] = {mpi_vec3_type, mpi_vec3_type, mpi_vec3_type, MPI_DOUBLE, MPI_INT, MPI_INT};
    MPI_Datatype mpi_body_type;
    MPI_Aint     offsets1[6];

    offsets1[0] = offsetof(body, position);
    offsets1[1] = offsetof(body, velocity);
    offsets1[2] = offsetof(body, force);
    offsets1[3] = offsetof(body, mass);
    offsets1[4] = offsetof(body, type);
    offsets1[5] = offsetof(body, number);

    MPI_Type_create_struct(nitems1, blocklengths1, offsets1, types1, &mpi_body_type);
    MPI_Type_commit(&mpi_body_type);
    


    //initializing input
    Initialize_Input(argc, argv, &numParticlesLight, &numParticleMedium, &numParticleHeavy, &numSteps, &subSteps, &timeSubStep, &width, &height, &output);\
    //setting number of total particles
    total_particles = numParticlesLight+numParticleMedium+numParticleHeavy;
    //allocating memory for bodies
    bodies = (body*)malloc(total_particles*sizeof(body));
    //Initializing particles, so that we can calulate forces
    if(my_rank == MASTER){
        Initialize_particles(&bodies, total_particles, numParticlesLight, numParticleMedium, numParticleHeavy, width, height);

    }
    //rbuf for body that will be given to processes initially
    body* rbuf;
    rbuf = (body*)malloc((2*(total_particles/p))*sizeof(body));
    //cbuf for body that will be passed between processes
    body* cbuf;
    cbuf = (body*)malloc((2*(total_particles/p))*sizeof(body));
    
    //scattering rbuf and cbuf to processes
    MPI_Scatter(bodies, (total_particles/p), mpi_body_type, rbuf, (total_particles/p), mpi_body_type, MASTER, MPI_COMM_WORLD); 
    MPI_Scatter(bodies, (total_particles/p), mpi_body_type, cbuf, (total_particles/p), mpi_body_type, MASTER, MPI_COMM_WORLD); 
    
    //for number of steps
    for (currentstep = 0; currentstep < numSteps; currentstep++){
        if (my_rank == MASTER){
            //setting the output file's name so that it has format argv[9]dddd.bmp where dddd is number starting from 0000.
            char *prefix;
            if (filecounter < 10){
                prefix = (char* )"000";
            }
            else if ((filecounter >= 10) && (filecounter < 100)){
                prefix = (char* )"00";
            }
            else if ((filecounter >= 100) && (filecounter < 1000))
                prefix = (char* )"0";
            else
                prefix = (char* )"";
            asprintf(&filename,"%s%s%d",argv[9],prefix,filecounter);
            saveIt(&bodies, filename, width, height, total_particles);
        }
        //for number of substeps
        for (currentsubstep = 0; currentsubstep < subSteps; currentsubstep++){
            if (my_rank == MASTER)
                //start timer
                start = MPI_Wtime();
            for (i = 0; i < p; i++){
                if (p != 1){
                    if (my_rank != 0) {
                        //send cbuf to different process. process with number p sends it to proccess p+1
                        MPI_Send(cbuf, (total_particles/p), mpi_body_type, (my_rank + 1) % p, i, MPI_COMM_WORLD);
                        MPI_Recv(cbuf, (total_particles/p), mpi_body_type, my_rank - 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        
                    }
                

                    // Now process 0 can receive from the last process.
                    //Master process have mbuf which will recieve data to achive ring communication.
                    if (my_rank == 0) {
                        body* mbuf;
                        mbuf = (body*)malloc((2*(total_particles/p))*sizeof(body));
                        MPI_Recv(mbuf, (total_particles/p), mpi_body_type, p - 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(cbuf, (total_particles/p), mpi_body_type, (my_rank + 1) % p, i, MPI_COMM_WORLD);
                        cbuf = mbuf;
                    }

                }
                //calculate the force among all processes.
                calculateForce(&rbuf, &cbuf,(total_particles/p));
                    
                

            }
            //update the position
            updatePosVel(&rbuf, &cbuf, (total_particles/p), timeSubStep, width, height);
            //gather the rbuf from all processes to bodies to plot output file
            MPI_Gather(rbuf, (total_particles/p), mpi_body_type, bodies, (total_particles/p), mpi_body_type, MASTER, MPI_COMM_WORLD);
            if (my_rank == MASTER){
                end = MPI_Wtime();
                //for calculating shortest and latest time
                if ((end-start) < fastest)
                    fastest = (end-start);
                if ((end-start) > latest)
                    latest = (end-start);
                total_time += (end-start);
            }
            
        }
        //increment counter for file name
        if (my_rank == 0){
            filecounter += 1;
        }

    }
    //print for result
    if (my_rank == MASTER)
    printf("<%f> <%f> <%f> \n", latest, fastest, (total_time/(numSteps*subSteps)));

    //free the allocated stuff
    MPI_Type_free(&mpi_body_type);
    MPI_Type_free(&mpi_vec3_type);
    free(bodies);
    MPI_Finalize();
    return 0;
}

void Initialize_Input(int argc, char* argv[], int* numParticlesLight, int* numParticleMedium, int* numParticleHeavy, int* numSteps, int* subSteps, double* timeSubStep, int* width, int* height, char** output){
    if(my_rank==MASTER){
        if( argc != 10){
            printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
            MPI_Finalize();
            exit(0);
        }
        
        *numParticlesLight = strtol(argv[1], NULL, 10);
        *numParticleMedium = strtol(argv[2], NULL, 10);
        *numParticleHeavy = strtol(argv[3], NULL, 10);
        *numSteps = strtol(argv[4], NULL, 10);
        *subSteps = strtol(argv[5], NULL, 10);
        *timeSubStep = strtod(argv[6], NULL);
        *width = strtol(argv[7], NULL, 10);
        *height = strtol(argv[8], NULL, 10);
        
    }
    //If the number of particles are negative
    if (numParticlesLight<0 || numParticleMedium<0 || numParticleHeavy<0){
        if (my_rank == MASTER){
            printf("Number of particles can't be negative \n");
        }
        MPI_Finalize();
        exit(0);
    }
    //If the number of steps or number of subSteps is negative
    if (numSteps<=0 || subSteps<=0){
        if (my_rank == MASTER){
            printf("Number of steps must be greater than zero \n");
        }
        MPI_Finalize();
        exit(0);
    }
    //If Time of substeps is negative
    if (timeSubStep<=0 ){
        if (my_rank == MASTER){
            printf("timeSubStep must be greater than zero \n");
        }
        MPI_Finalize();
        exit(0);
    }
    //If width or height is negative
    if (width<=0 || height<=0){
        if (my_rank == MASTER){
            printf("Width and height must be greater than zero \n");
        }
        MPI_Finalize();
        exit(0);

    }
    
    

    //Bcase to other processes
    MPI_Bcast(numParticlesLight, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numParticleMedium, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numParticleHeavy, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numSteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(subSteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(timeSubStep, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(width, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(height, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    
}


void Initialize_particles(body **bodies, int total_particles, int numParticlesLight, int numParticleMedium, int numParticleHeavy, int width, int height){
    if (my_rank == MASTER){
        double temp, temp1;
        int counter = 0;
        int i;
        int direction;
        srand48(time(NULL));
        //Light particles
        for (i = 0; i < numParticlesLight; i++){
            (*bodies)[i].type = 0;
            temp = drand48();
            temp = drand48();
            (*bodies)[i].position.x = (temp*width);
            temp = drand48();
            (*bodies)[i].position.y = (temp*height);
            (*bodies)[i].position.z = 0;
            temp = drand48();
            temp1 = drand48();
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.x = (direction*(temp*(((velocityLightMax-velocityLightMin)+velocityLightMin)/sqrt(2.0))));
            temp = drand48();
            temp1 = drand48();
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.y = (direction*(temp*(((velocityLightMax-velocityLightMin)+velocityLightMin)/sqrt(2.0))));
            (*bodies)[i].velocity.z = 0;
            temp = drand48();
            (*bodies)[i].mass = ((temp*(massLightMax-massLightMin)+massLightMin));
            (*bodies)[i].force.x = 0;
            (*bodies)[i].force.y = 0;
            (*bodies)[i].force.z = 0;
            (*bodies)[i].number = i;
        }
        //Medium particles
        counter = numParticlesLight;
        for (i = counter; i < numParticleMedium+counter; i++){
            (*bodies)[i].type = 1;
            temp = drand48();
            temp = drand48();
            (*bodies)[i].position.x = (temp*width);
            temp = drand48();
            (*bodies)[i].position.y = (temp*height);
            (*bodies)[i].position.z = 0;
            temp = drand48();
            temp1 = drand48();
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.x = (direction*(temp*(((velocityMediumMax-velocityMediumMin)+velocityMediumMin)/sqrt(2.0))));
            temp = drand48();
            temp1 = drand48();
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.y = (direction*(temp*(((velocityMediumMax-velocityMediumMin)+velocityMediumMin)/sqrt(2.0))));
            (*bodies)[i].velocity.z = 0;
            temp = drand48();
            (*bodies)[i].mass = ((temp*(massMediumMax-massMediumMin)+massMediumMin));
            (*bodies)[i].force.x = 0;
            (*bodies)[i].force.y = 0;
            (*bodies)[i].force.z = 0;
            (*bodies)[i].number = i;
        }
        //Heavy particles
        counter += numParticleMedium;
        for (i = counter; i < numParticleHeavy+counter; i++){
            (*bodies)[i].type = 2;
            temp = drand48();
            temp = drand48();
            (*bodies)[i].position.x = (temp*width);
            temp = drand48();
            (*bodies)[i].position.y = (temp*height);
            (*bodies)[i].position.z = 0;
            temp = drand48();
            temp1 = drand48();
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.x = (direction*(temp*(((velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin)/sqrt(2.0))));
            temp = drand48();
            temp1 = drand48();
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.y = (direction*(temp*(((velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin)/sqrt(2.0))));
            (*bodies)[i].velocity.z = 0;
            temp = drand48();
            (*bodies)[i].mass = ((temp*(massHeavyMax-massHeavyMin)+massHeavyMin));
            (*bodies)[i].force.x = 0;
            (*bodies)[i].force.y = 0;
            (*bodies)[i].force.z = 0;
            (*bodies)[i].number = i;
        }
        counter += numParticleHeavy;

    }
    
}
//for debugging
void print_status(body **bodies, int i){
    printf("particle %d has type : %d with x position : %f, y position : %f, x velocity : %f, y velocity : %f, mass : %f, xforce : %.20f, yforce : %.20f \n", i, (*bodies)[i].type, (*bodies)[i].position.x, (*bodies)[i].position.y, (*bodies)[i].velocity.x, (*bodies)[i].velocity.y, (*bodies)[i].mass, (*bodies)[i].force.x, (*bodies)[i].force.y);
    
}
//saving the image file
void saveIt(body** bodies, char *filename, int width, int height, int total_particles){
    int i;
    unsigned char* image;
    int holder = height*width;
    int tempx, tempy, temptype;

    image =(unsigned char*) malloc(3*holder); 

    //set the background to black
    if (my_rank == MASTER){
        for (i=0;i<holder;++i){
            image[3*i+0]=0;
            image[3*i+1]=0;
            image[3*i+2]=0;
        }
        //plot particles 
        for (i = 0; i < total_particles; i++){
            tempx = (*bodies)[i].position.x;
            tempy = (*bodies)[i].position.y;
            temptype = (*bodies)[i].type;
            vec3 tempcolor;
            if (temptype == 0)
                tempcolor = colourLight;
            else if (temptype==1)
                tempcolor = colourMedium;
            else
                tempcolor = colourHeavy;
            //if the particle is out of the box, make it invisible so that it won't cause error
            if (tempx >= width){
                tempcolor = colourNothing;
                tempx = width-1;
            }
            if (tempx <= 0){
                tempcolor = colourNothing;
                tempx = 0;
            }
            if (tempy >= height){
                tempcolor = colourNothing;
                tempy = height-1;
            }
            if (tempy <= 0){
                tempcolor = colourNothing;
                tempy = 0;
            }
            
            image[3*(tempy*width+tempx)+0]=(int(tempcolor.x)*255);
            image[3*(tempy*width+tempx)+1]=(int(tempcolor.y)*255);
            image[3*(tempy*width+tempx)+2]=(int(tempcolor.z)*255);
            
        
        }
        strcat(filename,".bmp");
        saveBMP(filename, image, width, height);
    }
    free(image);
}

void calculateForce(body** body1, body** body2, int a){
    int first, second;
    double temp, tempx, tempy;
    double distance, distancex;
    double angle;
    //if number of particle from body2 is greater than particle from body1, then calculate the force
    for (first = 0; first < a; first++){
        for (second = 0; second < a; second++){
            if ((*body2)[second].number > (*body1)[first].number){
                distance = sqrt(pow(fabs((*body1)[first].position.x - (*body2)[second].position.x), 2.0)+pow(fabs((*body1)[first].position.y - (*body2)[second].position.y), 2.0));
                distancex = fabs((*body1)[first].position.x - (*body2)[second].position.x);
                angle = acos((distancex/distance));
                temp = (G*((*body1)[first].mass)*((*body2)[second].mass))/(pow(distance, 2.0));
                tempx = (temp*cos(angle));
                tempy = (temp*sin(angle));
                if ((*body1)[first].position.x <= (*body2)[second].position.x){
                    (*body1)[first].force.x += tempx;
                    (*body2)[second].force.x += (-1 * tempx);
                }
                else{
                    (*body1)[first].force.x += (-1 * tempx);
                    (*body2)[second].force.x += tempx;
                }
                
                tempy = (G*((*body1)[first].mass)*((*body2)[second].mass))/(pow(distance,2.0));
                if ((*body1)[first].position.y <= (*body2)[second].position.y){
                    (*body1)[first].force.y += tempy;
                    (*body2)[second].force.y += (-1 * tempy);
                }
                else{
                    (*body1)[first].force.y += (-1 * tempy);
                    (*body2)[second].force.y += tempy;
                }
            }
        }
    }
   
}

//update the position and the velocity of the particles
void updatePosVel(body** bodies,body** body1, int total_particles, double timeSubStep, int width, int height){
    int a;
    int thistype;
    double xacceleration, yacceleration;
    double max, min;
    for (a = 0; a < total_particles; a++){
        (*bodies)[a].force += (*body1)[a].force;
        xacceleration = ((*bodies)[a].force.x)/((*bodies)[a].mass);
        yacceleration = ((*bodies)[a].force.y)/((*bodies)[a].mass);
        (*bodies)[a].velocity.x += (xacceleration*timeSubStep);
        (*bodies)[a].velocity.y += (yacceleration*timeSubStep);
        (*bodies)[a].velocity.z = 0;
        thistype = (*bodies)[a].type;
        
        if (thistype == 0){
            max = velocityLightMax;
            min = velocityLightMin;
        }
        else if (thistype == 1){
            max = velocityMediumMax;
            min = velocityMediumMin;
        }
        else{
            max = velocityHeavyMax;
            min = velocityHeavyMin;
        }
        //if the particle is slower thatn minimum velocity, set it to minimum velocity
        if ((*bodies)[a].velocity.Magnitude()<min){
            (*bodies)[a].velocity.Normalize();
            (*bodies)[a].velocity.x*=min;
            (*bodies)[a].velocity.y*=min;
            (*bodies)[a].velocity.z = 0;
        }
        //if it is faster
        else if ((*bodies)[a].velocity.Magnitude()>max){
            (*bodies)[a].velocity.Normalize();
            (*bodies)[a].velocity.x*=max;
            (*bodies)[a].velocity.y*=max;
            (*bodies)[a].velocity.z = 0;
        }

        vec3 temptemp;
        temptemp.x = (*bodies)[a].position.x + ((*bodies)[a].velocity.x*timeSubStep);
        temptemp.y = (*bodies)[a].position.y + ((*bodies)[a].velocity.y*timeSubStep);
        temptemp.z = 0;
        (*bodies)[a].position = temptemp;
        
        //set forces to 0
        (*bodies)[a].force.x=0;
        (*bodies)[a].force.y=0;
        (*bodies)[a].force.z=0;

    }
}

