/* 


  Mesh Structures that will be used in the project 
  
  
  Provided functions : 
       
  

*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <libmesh6.h>
#include <lplib3.h>
#include <math.h>



typedef double double4d[4];
typedef double double3d[3];
typedef double double2d[2];
typedef int    int5d[5];
typedef int    int4d[4]; 
typedef int    int3d[3];
typedef int    int2d[2];
typedef int    int1d;
typedef unsigned long long int u64;
typedef long long int lint1d;


/* 
  
  A simple mesh structure to hold the geomety/connectivity 

  Provided function msh_init, msh_read, msh_write, msh_check
  The mesh files are based on the meshb format (https://github.com/LoicMarechal/libMeshb), can be viewed with ViZiR4, see README to install it. 

  The following functions has to be implemented  : 
    msh_boundingbox : compute the bounding box of mesh
    msh_neighbors   : build the set of triangles surrounding elements, efficient implementation with Hash Tab
    msh_neighborsQ2 : a quadratic setting of the neigbhoring structures 

*/


typedef struct t_mesh
{
  int Dim;
  int NbrVer, NbrTri, NbrEfr, NbrEdg,NbrVerMax,NbrTriMax;
  
  double4d Box;           // domain bounding box
  
  //--- Data for the list of vertices
  double2d *Crd;          // coordinates
  
  //--- Data for the list of triangles
  int3d    *Tri;          // indices of the three vertices composing the triangle
  int3d    *TriVoi;       // indices of the three triangle neighbors
  int1d    *TriRef;       // triangle reference
  
  //--- Data for the list of boundary edges  
  int2d    *Efr;          // indices of the two vertices composing the boundary edge
  int2d    *EfrVoi;       // indices of the two boundary edges neighbors
  int1d    *EfrRef;       // boundary edges reference

  //--- Data for the list of edges  
  int2d    *Edg;          // indices of the two vertices composing the edge
    
} Mesh; 



//--- Provided functions 
Mesh   * msh_init();
Mesh   * msh_read(char *file, int readEfr);
int      msh_write(Mesh *Msh, char *file); 
double * sol_read(char *file, int mshDim, int mshNbrSol);

double   area(Mesh *Msh, int iTri);
double   inradius(Mesh *Msh, int iTri);
double   q1_quality(Mesh *Msh, int iTri);
double   q2_quality(Mesh *Msh, int iTri);

typedef struct {
  double val1;
  double val2;
  double val3;
} troisValeurs;

typedef struct  {
  double val1;
  double val2;
} deuxValeurs;

typedef struct  {
  int val1;
  int val2;
} deuxEntiers;

typedef struct {
  int nbaretetotal;
  int nbfrontiere;
  int nbaretevoisin;
} NeighborResult;

typedef struct node
{
    int data;
    struct node *next;
} Node;

typedef struct nodeEdg
{
    deuxEntiers data;
    struct nodeEdg *next;
} NodeEdg;

//--- Functions to be implemented  
int    msh_boundingbox(Mesh *Msh);         // compute the bouding box of the mesh            
NeighborResult    msh_neighbors(Mesh *Msh);           // build TriVoi with a hash table                 
int    msh_neighborsQ2(Mesh *Msh);         // build TriVoi with the naive quadratic approach 
  

//int nb_collisions(HashTable* hsh);

//--- A provided simple hash table data structure 
typedef struct hash_table
{
  int   SizHead;      // Maximum key value, the key is in [0,SizHead-1]
  int   NbrObj;       // Number of objects in the hash table
  int   NbrMaxObj;    // Maximum number of objects that can be store in the hash tab 
  int   *Head ;       // Head[key%(SizHead)] = link to the first object having this key in the LstObj list 
  int5d *LstObj;      // List of objects in the hash table 
  
  // LstObj[id][0:1] = iVer1-iVer2, the two vertices defining the edge  
  // LstObj[id][2:3] = iTri1,iTri2, the two neighboring triangles having iVer1-iVer2 as vertices 
  // LstObj[id][4]   = idxNxt,       the link to the next element in collision, if = 0 last element of the list 
  
} HashTable;


//--- Implementing the following function should be necessary 
HashTable * hash_init(int SizHead, int NbrMaxObj);          // alloc and set htable ==> allocate Head, LstObj 
int hash_key(HashTable *hsh, int iVer1 , int iVer2);
int hash_find(HashTable *hsh, int iVer1, int iVer2);            // return the id found (in LstObj ), if 0 the object is not in the list 
int hash_add (HashTable *hsh, int iVer1, int iVer2, int iTri);  // ==> add this entry in the hash tab 
int hash_suppr(HashTable *hsh, int iVer1, int iVer2, int iTri);  // ==> suppress this entry in the hash tab 

//HashTable genere_table(Mesh *Msh);
int aretefrontiere(Mesh *Msh);
int collisions(Mesh *Msh);

//--- Fonction used for adaptation 


//--- Writes a 2d metric tensor field in a file (extension should be .sol)
int msh_write2dmetric(char *file, int nmetric, double3d *metric);  
int msh_write2dfield_Triangles(char *file, int nfield, double *field);
int msh_write2dfield_Vertices(char *file, int nfield, double *field);


/////////////////////////PARTIE2///////////////////////////////
int localiser(Mesh *msh, double* P);
double areasigned(double* P1, double* P2, double *P3);
deuxValeurs circumcentre(double* P1, double* P2, double* P3);
troisValeurs longueur(double* P1, double* P2, double* P3);
int VerifyInCircle(Mesh *msh, int iTri, double *P);
int empiler(Node **PILE, int val);
int empiler_edge(NodeEdg **PILE, deuxEntiers val);
int depiler(Node **PILE, int *val);
int depilerEdg(NodeEdg **PILE, deuxEntiers *val);
int InPile(Node **PILE, int val);
int InPileEdg(NodeEdg **PILE, deuxEntiers val);
Node *cavity(Mesh *msh, double *P);
void print_pile(Node *pile);
NodeEdg *frontiere(Mesh *msh, Node* cavité);
int delaunay(Mesh *Msh, double *P);
void print_pile_Edge(NodeEdg *pile);
int msh_split_triangle(Mesh *Msh, int iTri, double *P);
int efface_contenu(Mesh *Msh, double *P);

void detect_contours_sobel(Mesh* msh, double* sol, double seuil, double2d *ListeP, double* Sol_Com, int* n);
double PSNR(Mesh *msh_compressed, Mesh *msh_origine, double* sol_compressed, double* sol_originale);