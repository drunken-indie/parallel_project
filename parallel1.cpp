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
#define G 6.67408*10000000 // The gravitational constant

int my_rank;

void Initialize_Input(int argc, char* argv[], int* numParticlesLight, int* numParticleMedium, int* numParticleHeavy, int* numSteps, int* subSteps, double* timeSubStep, int* width, int* height, char** output);
void Initialize_particles(body** bodies, int total_particles, int numParticlesLight, int numParticleMedium, int numParticleHeavy, int width, int height);
void print_status(body** bodies, int i);
void saveIt(body** bodies, char *filename, int width, int height, int total_particles);
void calculateForce(body** bodies, int a, int b);
void updatePosVel(body** bodies, int total_particles, double timeSubStep, int width, int height);
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
    struct body *bodies;
    struct body *processbody;
    char* filename;
    int filecounter=0;
    double start, end;
    double total_time=0;
    int i, j, k;
    body recv;
    
    

    
    //width=atoi(argv[7]);
    //height=atoi(argv[8]);

    
    
    //printf("sdsd \n"); 
    
    
    
    MPI_Init(&argc,&argv);

    int p;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    

    //printf("%d \n", p);

    /*create a tyoe for struct vec3*/
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


    /*create a tyoe for struct vec3*/
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

    //struct mpi_body_type* bodies;
    //root node stuff goes here
    Initialize_Input(argc, argv, &numParticlesLight, &numParticleMedium, &numParticleHeavy, &numSteps, &subSteps, &timeSubStep, &width, &height, &output);
    total_particles = numParticlesLight+numParticleMedium+numParticleHeavy;
    MPI_Bcast(&total_particles, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    
    
    bodies = (body*)malloc(total_particles*sizeof(body));
    
        
    
    if (my_rank == MASTER){
        Initialize_particles(&bodies, total_particles, numParticlesLight, numParticleMedium, numParticleHeavy, width, height);
        MPI_Bcast(&bodies, total_particles, mpi_body_type, MASTER, MPI_COMM_WORLD);

    }
    


    /*

    if (my_rank == MASTER){
        int to = 0;
        int what = 0;
        for (i = 0; i <total_particles; i++){
            to = i%p;
            MPI_Send(&bodies[i], 1, mpi_body_type, to, 0, MPI_COMM_WORLD);
            printf("Rank %d: sent %d structure body to %d Rank\n", my_rank, i, to);
            //printf("%d \n", to);
        }
        
        
    }
    if (my_rank >= 1){
        printf("Not Master \n");
        int j;
        MPI_Status status;
        body recv;
        for (j = 0; j < (total_particles/p); j++){
            body recv;
            MPI_Recv(&recv, 1, mpi_body_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Rank %d: Received: number = %d type = %d\n", my_rank,
                 recv.number, recv.type);
        }
        
    }

*/


/*MASTER TASK*/
    /*
    if (my_rank == MASTER){
        int to = 0;
        int what = 0;
        for (i = 1; i < p; i++){
            loop = (total_particles/p);
            MPI_Send(&loop, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        for (i = 0; i <total_particles; i++){
            to = i%p;
            MPI_Send(&bodies[i], 1, mpi_body_type, to, 0, MPI_COMM_WORLD);
            printf("Rank %d: sent %d structure body to %d Rank\n", my_rank, i, to);
        }
        body recv;
        for (j = 0; j < loop; j++){
            MPI_Recv(&recv, 1, mpi_body_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Rank %d: Received: number = %d type = %d\n", my_rank,
                 recv.number, recv.type);
        }
        
    }

*/
/*WORKER TASK*/
    /*
    else{
        MPI_Recv(&loop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("%d \n", loop);
        for (k = 0; k < loop; k++){
            MPI_Recv(&recv, 1, mpi_body_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Rank %d: Received: number = %d type = %d\n", my_rank,
                 recv.number, recv.type);
        }
    }
    */

    //double currentTime = 0;
    //printf("numSteps : %d \n", numSteps);
    for (currentstep = 0; currentstep < numSteps; currentstep++){
        
        if (my_rank == MASTER){
            char *prefix;
            //printf("Setting name %d \n", filecounter);
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
            //printf("Saving \n");
            saveIt(&bodies, filename, width, height, total_particles);
        }
            
  		
        //printf("%d \n", image[100]);
        for (currentsubstep = 1; currentsubstep <= subSteps; currentsubstep++){
            if (my_rank == MASTER){
                //print_status(&bodies, 0);
                //print_status(&bodies, 1);
            }
            
            //currentTime = (currentstep*subSteps+currentsubstep)*timeSubStep;
            start = MPI_Wtime();
            //printf("Calculating forces \n");
            /*MASTER TASK*/
            for (i = 0; i < p; i++){
                if (i == 0){
                    if (my_rank == MASTER){
                        int to = 0;
                        for (i = 0; i <total_particles; i++){
                            to = i%p;
                            MPI_Send(&bodies[i], 1, mpi_body_type, to, 0, MPI_COMM_WORLD);
                            //printf("Rank %d: sent %d structure body to %d Rank\n", my_rank, i, to);
                        }
                    }


                    /*WORKER TASK*/
                    
                    if (my_rank >= MASTER){
                        processbody = (body*)malloc((2*(total_particles/p))*sizeof(body));
                        for (k = 0; k < (total_particles/p); k++){
                            MPI_Recv(&recv, 1, mpi_body_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            //printf("Rank %d: Received: number = %d type = %d\n", my_rank,
                              // recv.number, recv.type);
                            processbody[k] = recv;
                            //printf("recv : %d \n", recv.number);
                            
                        }
                        
                        //printf("a: %d \n", processbody[1].number);
                        //print_status(&processbody, 1);
                        MPI_Barrier(MPI_COMM_WORLD);
                        //CALCULATEEEE
                        calculateForce(&processbody, 0, 2);
                        //printf("a: %d \n", processbody[0].type);
                        if (my_rank == MASTER){
                           // print_status(&processbody, 0);
                           // print_status(&processbody, 1);
                        }
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    
            
            
                }
                else{
                        
                    if (my_rank == MASTER){
                        printf("ahah \n");
                    }
                }
            }
            


            //calculateForce(&bodies, total_particles);
            //printf("Updating pos vel \n");
            updatePosVel(&bodies, total_particles, timeSubStep, width, height);
            end = MPI_Wtime();
            total_time += (end-start);
        }
        filecounter += 1;
        //printf("filename : %s \n", filename);
            
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == MASTER){
        printf("Total Time : %f \n", total_time);
    }
    
        
    
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
    
    MPI_Bcast(numParticlesLight, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numParticleMedium, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numParticleHeavy, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numSteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(subSteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(timeSubStep, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(width, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(height, 1, MPI_CHAR, MASTER, MPI_COMM_WORLD);

    /*
    MPI_Bcast(numParticlesLight, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numParticleMedium, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numParticleHeavy, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(numSteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(subSteps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(timeSubStep, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(width, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(height, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    */
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
            //printf("temp1 : %d \n", int(temp1+0.5));
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.x = (direction*(temp*(velocityLightMax-velocityLightMin)+velocityLightMin));
            temp = drand48();
            temp1 = drand48();
            //printf("temp1 : %d \n", int(temp1+0.5));
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.y = (direction*(temp*(velocityLightMax-velocityLightMin)+velocityLightMin));
            (*bodies)[i].velocity.z = 0;
            temp = drand48();
            (*bodies)[i].mass = ((temp*(massLightMax-massLightMin)+massLightMin));
            (*bodies)[i].force.x = 0;
            (*bodies)[i].force.y = 0;
            (*bodies)[i].force.z = 0;
            (*bodies)[i].number = i;
            //printf("%f \n", temp);
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
            //printf("temp1 : %d \n", int(temp1+0.5));
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.x = (direction*(temp*(velocityMediumMax-velocityMediumMin)+velocityMediumMin));
            temp = drand48();
            temp1 = drand48();
            //printf("temp1 : %d \n", int(temp1+0.5));
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.y = (direction*(temp*(velocityMediumMax-velocityMediumMin)+velocityMediumMin));
            (*bodies)[i].velocity.z = 0;
            temp = drand48();
            (*bodies)[i].mass = ((temp*(massMediumMax-massMediumMin)+massMediumMin));
            (*bodies)[i].force.x = 0;
            (*bodies)[i].force.y = 0;
            (*bodies)[i].force.z = 0;
            (*bodies)[i].number = i;
            //printf("%f \n", temp);
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
            //printf("temp1 : %d \n", int(temp1+0.5));
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.x = (direction*(temp*(velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin));
            temp = drand48();
            temp1 = drand48();
            //printf("temp1 : %d \n", int(temp1+0.5));
            //positive
            if (int(temp1+0.5) == 0){
                direction = 1;
            }
            else
                direction = -1;
            (*bodies)[i].velocity.y = (direction*(temp*(velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin));
            (*bodies)[i].velocity.z = 0;
            temp = drand48();
            (*bodies)[i].mass = ((temp*(massHeavyMax-massHeavyMin)+massHeavyMin));
            (*bodies)[i].force.x = 0;
            (*bodies)[i].force.y = 0;
            (*bodies)[i].force.z = 0;
            (*bodies)[i].number = i;
            //printf("%f \n", temp);
        }
        counter += numParticleHeavy;

        //printf("%d \n", counter);
    }
    
}

void print_status(body **bodies, int i){
    printf("particle %d has type : %d with x position : %f, y position : %f, x velocity : %f, y velocity : %f, mass : %f, xforce : %f, yforce : %f \n", i, (*bodies)[i].type, (*bodies)[i].position.x, (*bodies)[i].position.y, (*bodies)[i].velocity.x, (*bodies)[i].velocity.y, (*bodies)[i].mass, (*bodies)[i].force.x, (*bodies)[i].force.y);
    
}

void saveIt(body** bodies, char *filename, int width, int height, int total_particles){
    int i;
    unsigned char* image;
    int holder = height*width;
    int tempx, tempy, temptype;

    image =(unsigned char*) malloc(3*holder); 

    if (my_rank == MASTER){
        for (i=0;i<holder;++i){
            image[3*i+0]=0;
            image[3*i+1]=0;
            image[3*i+2]=0;
        }
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

void calculateForce(body** bodies, int a, int b){
    int first, second;
    double temp, tempx, tempy;
    double distance, distancex;
    double angle;
    for (first = 0; first < b; first++){
        for (second = 0; second < b; second++){
            if ((*bodies)[second].number > (*bodies)[first].number){
                distance = sqrt(pow(fabs((*bodies)[first].position.x - (*bodies)[second].position.x), 2.0)+pow(fabs((*bodies)[first].position.y - (*bodies)[second].position.y), 2.0));
                distancex = fabs((*bodies)[first].position.x - (*bodies)[second].position.x);
                angle = acos((distancex/distance));
                //printf("angle : %f, aa : %f bb : %f \n", angle, cos(angle), (distancex/distance));
                temp = (G*((*bodies)[first].mass)*((*bodies)[second].mass))/(pow(distance, 2.0));
                tempx = (temp*cos(angle));
                tempy = (temp*sin(angle));
                //printf("xdistance : %f \n", distance);
                //printf("%f \n", (G*(bodies[first].mass)*(bodies[second].mass)));
                //printf("tempx : %f \n", tempx);
                if ((*bodies)[first].position.x <= (*bodies)[second].position.x){
                    (*bodies)[first].force.x += tempx;
                    (*bodies)[second].force.x += (-1 * tempx);
                }
                else{
                    (*bodies)[first].force.x += (-1 * tempx);
                    (*bodies)[second].force.x += tempx;
                }
                
                //printf("ydistance : %f \n", distance);
                //printf("%f \n", (G*(bodies[first].mass)*(bodies[second].mass)));
                tempy = (G*((*bodies)[first].mass)*((*bodies)[second].mass))/(pow(distance,2.0));
                //printf("tempy : %f \n", tempy);
                if ((*bodies)[first].position.y <= (*bodies)[second].position.y){
                    (*bodies)[first].force.y += tempy;
                    (*bodies)[second].force.y += (-1 * tempy);
                }
                else{
                    (*bodies)[first].force.y += (-1 * tempy);
                    (*bodies)[second].force.y += tempy;
                }
                
                //printf("first xforce : %f \t yforce : %f \n", bodies[first].force.x, bodies[first].force.y);
            }
        }
    }
}

void updatePosVel(body** bodies, int total_particles, double timeSubStep, int width, int height){
    int a;
    int thistype;
    double xacceleration, yacceleration;
    double max, min;
    for (a = 0; a < total_particles; a++){
        //print_status(bodies, a);
        xacceleration = ((*bodies)[a].force.x)/((*bodies)[a].mass);
        //printf("xacceleration : %f \t", xacceleration);
        yacceleration = ((*bodies)[a].force.y)/((*bodies)[a].mass);
        //printf("yacceleration : %f \n", yacceleration);
        (*bodies)[a].velocity.x += (xacceleration*timeSubStep);
        //printf("aa %f \t", (xacceleration*timeSubStep));
        (*bodies)[a].velocity.y += (yacceleration*timeSubStep);
        (*bodies)[a].velocity.z = 0;
        //printf("bb %f \n", (yacceleration*timeSubStep));
        thistype = (*bodies)[a].type;
        //printf("%d \n", thistype);
        
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
        //bodies[a].velocity.Normalize();
        //printf("%f \t", bodies[a].velocity.x);
        //printf("%f \+t", bodies[a].velocity.y);
        //printf("%f \n", bodies[a].velocity.z);
        
        if ((*bodies)[a].velocity.Magnitude()<min){
            //printf("a \n");
            (*bodies)[a].velocity.Normalize();
            //print_status(bodies, a);
            (*bodies)[a].velocity.x*=min;
            (*bodies)[a].velocity.y*=min;
            (*bodies)[a].velocity.z = 0;
        }
        
        else if ((*bodies)[a].velocity.Magnitude()>max){
            (*bodies)[a].velocity.Normalize();
            //print_status(bodies, a);
            (*bodies)[a].velocity.x*=max;
            (*bodies)[a].velocity.y*=max;
            (*bodies)[a].velocity.z = 0;
        }

        vec3 temptemp;
        temptemp.x = (*bodies)[a].position.x + ((*bodies)[a].velocity.x*timeSubStep);
        temptemp.y = (*bodies)[a].position.y + ((*bodies)[a].velocity.y*timeSubStep);
        temptemp.z = 0;
        (*bodies)[a].position = temptemp;
        
        //Border
        
        
        
        (*bodies)[a].force.x=0;
        (*bodies)[a].force.y=0;
        (*bodies)[a].force.z=0;

    }
}

