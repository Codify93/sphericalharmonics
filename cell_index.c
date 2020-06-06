#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spherical_harmonics.h"

#define N 12288
#define PI 3.14159265359

#define distance(x1,x2,y1,y2,z1,z2) sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2))


typedef struct{
    char a[2];
    double posx;
    double posy;
    double posz;
    int label;
    unsigned long int num;
    unsigned long int cell_index;
}particle;


typedef struct Node {
    particle member;
    struct Node * next_member;
}Node;

typedef struct {
    Node * cluster;
    unsigned long int m; // cell index
}Cell;

typedef struct{
    Node * head;
}Head;

Node * create_node();
Node * add_node(Node * head,Node * tail);

int main(){
    FILE * fp;
    FILE * fpw; 
    FILE * test;
    fp = fopen("asilica.txt", "r");
    fpw = fopen("output.txt", "w+");
    test = fopen("testing.txt", "w+");
    particle * mono;

    Head * heads;

    Node * temp_pointer;

    double theta;

    double phi;

    double mod_r;

    double rc = 4.0;
    double L = 60.0;
    int Lx, Ly, Lz;
    unsigned long int Cx, Cy, Cz, C, Lxyz, nL=0;

    int max = 0;

    Node * tpointer;


    Lx = Ly = Lz = floor(L/rc);

    Lxyz = (Lx+1)*(Ly+1)*(Lz+1);


    int totaln[Lxyz];


    for(unsigned int Ctemp = 0; Ctemp < Lxyz; Ctemp++){
            totaln[Ctemp] = 0;
    }
     Cell * cell;
     cell = (Cell *)malloc((Lxyz+1)*sizeof(Cell));


    mono = (particle *)malloc(N*sizeof(particle));

    heads = (Head *)malloc((Lxyz+1)*sizeof(Head));

    for(unsigned int Ctemp = 0; Ctemp < Lxyz; Ctemp++){
            (cell+Ctemp)->m = -1;
    }

    for(unsigned int i = 0; i < N; i++){
        fscanf(fp,"%c%c   %lf   %lf   %lf  %d     %lu\n",&(mono+i)->a[0],&(mono+i)->a[1],&(mono+i)->posx,
                                                    &(mono+i)->posy,&(mono+i)->posz,&(mono+i)->label,&(mono+i)->num);
    
    (mono+i)->posx = (mono+i)->posx + (double)30;
    (mono+i)->posy = (mono+i)->posy + (double)30;
    (mono+i)->posz = (mono+i)->posz + (double)30; 


    Cx = ceil((mono+i)->posx / rc) ;
    Cy = ceil((mono+i)->posy / rc) ;
    Cz = ceil((mono+i)->posz / rc) ;

    C = Cx*Ly*Lz + Cy*Lz + Cz;

    (mono+i)->cell_index = C;

    temp_pointer = (Node *)malloc(sizeof(Node));


    (cell+C)->cluster = create_node();

        (cell+C)->cluster->member.a[0] = (mono+i)->a[0];
        (cell+C)->cluster->member.a[1] = (mono+i)->a[1];
        (cell+C)->cluster->member.posx = (mono+i)->posx;
        (cell+C)->cluster->member.posy = (mono+i)->posy;
        (cell+C)->cluster->member.posz = (mono+i)->posz;
        (cell+C)->cluster->member.num  = (mono+i)->num;
        
        if((cell+C)->m == -1){
        (heads+C)->head = (cell+C)->cluster;
   
        (cell+C)->m = C;
        }
        else{
       
        (heads+C)->head = add_node((heads+C)->head, (cell+C)->cluster);
        totaln[C]++;
        }
    
    }
         
        unsigned long int C1;
        int cx,cy,cz;
        int tempC;

        double sum;
        double sq_average;
        int totalnn;
        int l = 4;
        double average;
       
       for(unsigned long int k = 0; k < N; k++){
           fprintf(test,"start"); 
       


            sum = 0;
            sq_average = 0;
            totalnn = 0;

            unsigned long int C1 = (mono+k)->cell_index ;

            cx = C1 / (Ly * Lz);
            cy = (C1 / Lz) % Ly;
            cz = C1 % Lz;

    

            for(int ml = -l; ml <= l; ml++){
               sum = 0;
               totalnn = 0;
            
                 for(int mx = -1; mx < 2; mx++)
                 for(int my = -1; my < 2; my++)
                 for(int mz = -1; mz < 2; mz++){
                  tempC = (cx+mx)*Ly*Lz + (cy+my)*Lz + (cz+mz);

                 fprintf(test,"tempC = %d\n",tempC);

                if(0 <= tempC && tempC < Lxyz+1){
                         tpointer = (heads+tempC)->head;
                    while(tpointer != NULL ){
                       
                mod_r = distance(tpointer->member.posx,(mono+k)->posx,tpointer->member.posy,
                (mono+k)->posy,tpointer->member.posz,(mono+k)->posz);

                   if(mod_r > 0 && mod_r < 2){



                    theta = atan((tpointer->member.posy - (mono+k)->posy)/ (tpointer->member.posx - (mono+k)->posx));
                     phi  = acos((tpointer->member.posz - (mono+k)->posz)/ mod_r);
                     sum = y_lm(phi, theta, l, ml) + sum ;

                     totalnn++ ; 
                
                   }

                tpointer = tpointer->next_member;
                    }
                }
            }   
                average = (double)sum / totalnn;
                sq_average = average*average + sq_average;
 
                
                }

                
         
        
          double Q_l = sqrt((4*PI*sq_average)/(2*l + 1));
         
                           fprintf(fpw,"%c%c   %lf    %lf     %lf     %lf              %lu     %lu\n",(mono+k)->a[0],
                           (mono+k)->a[1],(mono+k)->posx,(mono+k)->posy,(mono+k)->posz, Q_l, (mono+k)->num, (mono+k)->cell_index);
        
                            
       } 



       free(cell);
    free(mono);
    free(heads);

    return 0;
}



Node * create_node(){
    Node * temp;
    temp = (Node *)malloc(sizeof(Node));
    temp->next_member = NULL;
    return temp;
}

Node * add_node(Node * head,Node * tail){
    tail->next_member = head;
    return tail;
}



