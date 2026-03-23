#include <mesh.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[])
{ 
  int    iTri, iVer;
  double to, ti;
  bool hachage = false;
  char method[30] = "CompresserImageSobel";  //  "DelaunayAleatoire" ou "SplitSimple" ou "PetitTest" ou  "CompresserImageAleatoire" ou "CompresserImageGradLocal" ou "PSNR"
/*
PetitTest: test manuelle pour verifier le fonctionnement de la cavité et sa frontiere..
SplitSimple: pour diviser un triangle sur trois 
DelaunayAleatoire: pour ajouter des points aleatoirement en utilisant la cavité, Delaunay...
CompresserImageAleatoire: on utilise delaunay pour inserer des points 1 pixel/n pixels..
CompresserImageGradLocal: Compresser l'image en utilisant le gradient local...
CompresserImageSobel: Compresser l'image en utilisant Sobel...

*/


  if ( argc < 2 ) {
    printf(" usage : mesh file \n");
    return 0;
    }

  //--- read a mesh 
  to =  GetWallClock();
  Mesh * Msh = msh_read(argv[1], 0);
  ti =  GetWallClock();
  
  if ( ! Msh ) return 0;
  
FILE *fileTRInt = fopen("triinit.txt", "a"); // Ouvrir en mode ajout
for (int iEdg=0; iEdg<3; iEdg++){
  for (iTri=1; iTri<=Msh->NbrTri; iTri++){
    fprintf(fileTRInt,"  %d ;  ", Msh->Tri[iTri][iEdg]);
  }
  fprintf(fileTRInt, " \n");
}
fclose(fileTRInt);
  
  printf("  Vertices   %10d \n", Msh->NbrVer);
  printf("  Triangles  %10d \n", Msh->NbrTri);
  printf("  time to read the mesh %lg (s) \n",ti-to);
  
  /*
  //--- create neigbhors Q2 version 
  to =  GetWallClock();
  msh_neighborsQ2(Msh);
  ti =  GetWallClock();
  printf("  time q2 neigh.        %lg (s) \n",ti-to);
  */
  if (hachage== true){
  int nfrontiere = aretefrontiere(Msh);
  printf("  nb d'arêtes frontiere  : %d \n", nfrontiere);
  int ncollision=  collisions(Msh);
  printf("  nb de collisions  : %d \n", ncollision);


  //--- create neigbhors with hash table 
  to =  GetWallClock();
  NeighborResult result= msh_neighbors(Msh);
  ti =  GetWallClock();
  printf("  time hash tab neigh.  %lg (s) \n",ti-to);
  printf("  nb d'arêtes total  : %d \n", result.nbaretetotal);
  printf("  nb d'arêtes frontiere  : %d \n", result.nbfrontiere);
  printf("  nb d'arêtes voisins  : %d \n", result.nbaretevoisin);
    
  FILE *file1 = fopen("quality1.txt", "a"); // Ouvrir en mode ajout
  FILE *file2 = fopen("quality2.txt", "a"); // Ouvrir en mode ajout




  //--- TODO: compute mesh quality
  double  *Qal1 = (double  *)malloc(sizeof(double ) * (Msh->NbrTri+1));
  double  *Qal2 = (double  *)malloc(sizeof(double ) * (Msh->NbrTri+1));
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
  	 Qal1[iTri] = q1_quality(Msh, iTri);
     Qal2[iTri] = q2_quality(Msh, iTri);
     fprintf(file1," %lg \n",Qal1[iTri] );
     fprintf(file2," %lg \n",Qal2[iTri] );
  } 
  fclose(file1);
  fclose(file2);
 

  msh_write2dfield_Triangles("quality1.solb", Msh->NbrTri, Qal1);
  msh_write2dfield_Triangles("quality2.solb", Msh->NbrTri, Qal2);
  }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                          PARTIE 2
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
FILE *fileTRI = fopen("tri.txt", "a"); 
for (int iEdg=0; iEdg<3; iEdg++){
  for (iTri=1; iTri<=Msh->NbrTri; iTri++){
    fprintf(fileTRI,"%d ;  ", Msh->Tri[iTri][iEdg]);
  }
  fprintf(fileTRI, " \n");
}
fclose(fileTRI);
*/
msh_neighborsQ2(Msh);

FILE *fileVoisin = fopen("voisin.txt", "a");
for (int iEdg=0; iEdg<3; iEdg++){
  fprintf(fileVoisin, " iEdg = %d \n", iEdg);
  for (iTri=1; iTri<=Msh->NbrTri; iTri++){
    fprintf(fileVoisin,"  %d ,  ", Msh->TriVoi[iTri][iEdg]);
  }
  fprintf(fileVoisin, " \n");
}
fclose(fileVoisin);

/*
NeighborResult result= msh_neighbors(Msh);
FILE *fileVoisinHACH = fopen("voisinHACH.txt", "a"); // Ouvrir en mode ajout
for (int iEdg=0; iEdg<3; iEdg++){
  fprintf(fileVoisinHACH, " iEdg = %d \n", iEdg);
  for (iTri=1; iTri<=Msh->NbrTri; iTri++){
    fprintf(fileVoisinHACH,"  %d ,  ", Msh->TriVoi[iTri][iEdg]);
  }
  fprintf(fileVoisinHACH, " \n");
}
fclose(fileVoisinHACH);
*/

msh_boundingbox(Msh);
/*
printf("xmin = %f \n",Msh->Box[0]);
printf("xmax= %f \n",Msh->Box[1]);
printf("ymin= %f \n",Msh->Box[2]);
printf("ymax= %f \n" ,Msh->Box[3]);

printf("ETAT INITIAL");
printf("Triangles: \n");
for(int kTri=1; kTri <= Msh->NbrTri; kTri++){
  printf("triangle num %d : \n", kTri);
  printf("Vertices: %d ; %d ; %d  \n", Msh->Tri[kTri][0],Msh->Tri[kTri][1],Msh->Tri[kTri][2]);
}
*/

///////SPLIT LE TRIANGLE EN TROIS : FONCTIONNE ////

if (strcmp(method, "SplitSimple") == 0) {
  srand(7);
  int N=5;
  double xmin = Msh->Box[0], xmax = Msh->Box[1];  // Domaine en x
  double ymin = Msh->Box[2], ymax = Msh->Box[3];  // Domaine en y

  for (int i = 1; i <= N; i++) {
    double2d p;
    p[0] = (xmin+0.1) + (xmax - xmin) * ((double)rand() / RAND_MAX);
    p[1] = (ymin+0.1) + (ymax - ymin) * ((double)rand() / RAND_MAX);

    printf("Point %d : (%.6f, %.6f)\n", i , p[0], p[1]);
  
    int loc=localiser(Msh,p);
    msh_split_triangle(Msh, loc, p);
    printf("Triangles: \n");
    for(int kTri=1; kTri <= Msh->NbrTri; kTri++){
      printf("triangle num %d : \n", kTri);
      printf("Vertices: %d ; %d ; %d  \n", Msh->Tri[kTri][0],Msh->Tri[kTri][1],Msh->Tri[kTri][2]);
    }
  }  
  msh_write(Msh,"triangledivise.mesh");
}

//////UTILISATION DELAUNEY/////
if (strcmp(method, "PetitTest") == 0) {
  double2d p={0.04,0.549};
  int loc=localiser(Msh,p);
  printf("loc= %d \n", loc);
  printf("p[0]= %f \n",p[0]);
  printf("p[1]= %f \n",p[1]);
  printf("Voisins: %d , %d , %d \n", Msh->TriVoi[1244][0], Msh->TriVoi[1244][1], Msh->TriVoi[1244][2]);
  Node *cavité=cavity(Msh,p);
  printf("cavité : \n");
  print_pile(cavité);
  printf("frontiere de cavité manuel : \n");
  printf("frontiere de 4825 : %d,%d ,%d\n", Msh->TriVoi[4825][0],Msh->TriVoi[4825][1],Msh->TriVoi[4825][2]);
  printf("frontiere de 3865 : %d,%d ,%d\n", Msh->TriVoi[3865][0],Msh->TriVoi[3865][1],Msh->TriVoi[3865][2]);
  printf("frontiere de 4354 : %d,%d ,%d\n", Msh->TriVoi[4354][0],Msh->TriVoi[4354][1],Msh->TriVoi[4354][2]);
  printf("iTri= 4354 i= %d, j= %d, k=%d  \n",Msh->Tri[4354][0],Msh->Tri[4354][1],Msh->Tri[4354][2]);
  printf("iTri= 4825 i= %d, j= %d, k=%d \n",Msh->Tri[4825][0],Msh->Tri[4825][1],Msh->Tri[4825][2]);
  NodeEdg *edgCavity= frontiere(Msh,cavité);
  print_pile_Edge(edgCavity);
  delaunay(Msh,p);
  
  msh_write(Msh,"newmesh.mesh");
  printf("Triangles: \n");
  for(int kTri=1; kTri <= Msh->NbrTri; kTri++){
    printf("triangle num %d : \n", kTri);
    printf("Vertices: %d ; %d ; %d  \n", Msh->Tri[kTri][0],Msh->Tri[kTri][1],Msh->Tri[kTri][2]);
  }
  printf("\n Sommets: \n");
  for(int kVer=1; kVer <= Msh->NbrVer; kVer++){
    printf("Sommet %d : \n", kVer);
    printf("Coordonnees: %f ; %f  \n", Msh->Crd[kVer][0],Msh->Crd[kVer][1]);
  }
  
  
  }
  
if (strcmp(method,"DelaunayAleatoire") == 0) {
  srand(7);
  int N=atoi(argv[2]);
  printf(N);
  double xmin = Msh->Box[0], xmax = Msh->Box[1];  // Domaine en x
  double ymin = Msh->Box[2], ymax = Msh->Box[3];  // Domaine en y
  to =  GetWallClock();
  int TriAjoutée =0;

  for (int i = 1; i <= N; i++) 
  {
    double2d p;
    p[0] = (xmin+0.1)+ (xmax - xmin) * ((double)rand() / RAND_MAX);
    p[1] = (ymin+0.1) + (ymax - ymin) * ((double)rand() / RAND_MAX);


    //printf("Point %d : (%.6f, %.6f)\n", i , p[0], p[1]);
  
    int n= delaunay(Msh,p);
    
    TriAjoutée += n;

    
    printf("MISE A JOUR \n");
    printf("Triangles: \n");
    for(int kTri=1; kTri <= Msh->NbrTri; kTri++){
      printf("triangle num %d : \n", kTri);
      printf("Vertices: %d ; %d ; %d  \n", Msh->Tri[kTri][0],Msh->Tri[kTri][1],Msh->Tri[kTri][2]);
    }
    printf("\n Sommets: \n");
    for(int kVer=1; kVer <= Msh->NbrVer; kVer++){
      printf("Sommet %d : \n", kVer);
      printf("Coordonnees: %f ; %f  \n", Msh->Crd[kVer][0],Msh->Crd[kVer][1]);
    }
    }
  ti =  GetWallClock();
  printf("N= %d time delaunay.  %lg (s) \n", N, ti-to);
  //printf("Ajout de %d triangles \n", TriAjoutée);
  msh_write(Msh,"nouveaumaillage.mesh");

}

if (strcmp(method,"CompresserImageAleatoire")==0){
  Mesh * msh_joc = msh_read("data/joconde.mesh", 0);
  //msh_neighbors(msh_joc);
  double* sol_joc = sol_read("data/joconde.sol",msh_joc->Dim,msh_joc->NbrVer);
  //int TriAjoutée =0;
  printf("nbreVer = %d \n",msh_joc->NbrVer);
  int n=1;
  double  *Sol_Com = (double *)malloc(sizeof(double ) * (msh_joc->NbrVer +1));
  double2d *ListeP =calloc( (msh_joc->NbrVer +1), sizeof(double2d) );
  for (iVer=1; iVer<=msh_joc->NbrVer; iVer++){
    if (iVer % 10 == 1) //choix d'un pixel sur 5
    {
      ListeP[n][0]=msh_joc->Crd[iVer][0];
      ListeP[n][1]=msh_joc->Crd[iVer][1];
      Sol_Com[n+4]=sol_joc[iVer];
      printf("p= %f , %f \n ",ListeP[n][0],ListeP[n][1] );
      n++;
    }
  }
  int TriAjoutée=0;
  for (int i = 1; i <= n; i++) 
  {
    double2d p;
    p[0]=ListeP[i][0];
    p[1]=ListeP[i][1];
    int ndelaunay= delaunay(Msh,p);
    TriAjoutée += ndelaunay;
  }
  printf("TriAjoutée= %d \n",TriAjoutée);
  msh_write(Msh,"maillageCompressed.mesh");
  msh_write2dfield_Vertices("maillageCompressed.sol", Msh->NbrVer ,Sol_Com);  
}

if (strcmp(method,"CompresserImageGradLocal")==0){
  Mesh * msh_joc = msh_read("data/joconde.mesh", 0);
  //msh_neighbors(msh_joc);
  double* sol_joc = sol_read("data/joconde.sol",msh_joc->Dim,msh_joc->NbrVer);
  printf("(i,j) correspond de iVer=290 = %f , %f \n",msh_joc->Crd[290][0],msh_joc->Crd[290][1]);
  //int TriAjoutée =0;
  printf("nbreVer = %d \n",msh_joc->NbrVer);
  int n=1;
  double  *Sol_Com = (double *)malloc(sizeof(double ) * (msh_joc->NbrVer +1));
  double2d *ListeP =calloc( (msh_joc->NbrVer +1), sizeof(double2d) );

  double seuil = 10.0; 
  for (int i=1;i<=4;i++){Sol_Com[i] = sol_joc[i];}
for (int iVer = 5; iVer < msh_joc->NbrVer ; iVer++) {
    double dx = msh_joc->Crd[iVer+1][0] - msh_joc->Crd[iVer][0];
    double dy = msh_joc->Crd[iVer+1][1] - msh_joc->Crd[iVer][1];
    double dsol = sol_joc[iVer+1] - sol_joc[iVer];
    double grad = fabs(dsol) / sqrt(dx*dx + dy*dy); 

    if (grad > seuil) {
        ListeP[n][0] = msh_joc->Crd[iVer][0];
        ListeP[n][1] = msh_joc->Crd[iVer][1];
        //Sol_Com[n+4] = sol_joc[iVer];
        Sol_Com[n] = sol_joc[iVer];
        n++;
    }
}
  printf("n= %d \n",n);
  int TriAjoutée =0;
  for (int i = 1; i <= n; i++) 
  {
    double2d p;
    p[0]=ListeP[i][0];
    p[1]=ListeP[i][1];
    int ndelaunay= delaunay(Msh,p);
    TriAjoutée += ndelaunay;
  }
  printf("TriAjoutée= %d \n",TriAjoutée);
  msh_write(Msh,"maillageCompressedGrad.mesh");
  msh_write2dfield_Vertices("maillageCompressedGrad.sol", Msh->NbrVer ,Sol_Com);
  //Verifier PSNR
double psnr= PSNR(Msh,msh_joc,Sol_Com,sol_joc);
//printf("psnr = %f \n ",psnr);
double2d p={1.0 , 440.0};
printf("loc= %d \n", localiser(Msh,p));

} 
FILE *fileTRI = fopen("tri.txt", "a"); 
for (int iEdg=0; iEdg<3; iEdg++){
  for (iTri=1; iTri<=43; iTri++){
    fprintf(fileTRI,"%d ;  ", Msh->Tri[iTri][iEdg]);
  }
  fprintf(fileTRI, " \n");
}
fclose(fileTRI);


if (strcmp(method,"CompresserImageSobel")==0){
  Mesh * msh_joc = msh_read("data/joconde.mesh", 0);
  //msh_neighbors(msh_joc);
  double* sol_joc = sol_read("data/joconde.sol",msh_joc->Dim,msh_joc->NbrVer);
  //int TriAjoutée =0;
  printf("nbreVer = %d \n",msh_joc->NbrVer);
  int n=1;
  double  *Sol_Com = (double *)malloc(sizeof(double ) * (msh_joc->NbrVer +1));
  double2d *ListeP =calloc( (msh_joc->NbrVer +1), sizeof(double2d) );
  double seuil = 40.0; 
  detect_contours_sobel(msh_joc, sol_joc, seuil, ListeP, Sol_Com, &n);
  printf("n= %d \n",n);
  int TriAjoutée =0;
  for (int i = 1; i <= n; i++) 
  {
    double2d p;
    p[0]=ListeP[i][0];
    p[1]=ListeP[i][1];
    int ndelaunay= delaunay(Msh,p);
    TriAjoutée += ndelaunay;
  }
  printf("TriAjoutée= %d \n",TriAjoutée);
  msh_write(Msh,"maillageCompressedSobel.mesh");
  msh_write2dfield_Vertices("maillageCompressedSobel.sol", Msh->NbrVer ,Sol_Com);
  
}







//--- TODO: compute metric field
  double3d *Met = (double3d *)malloc(sizeof(double3d) * (Msh->NbrVer+1));

  for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
  	 Met[iVer][0] = 1.;
  	 Met[iVer][1] = 0.;
  	 Met[iVer][2] = 1.;
  } 
  
  msh_write2dmetric("metric.solb", Msh->NbrVer, Met);
  	
  	
  return 0;
}






