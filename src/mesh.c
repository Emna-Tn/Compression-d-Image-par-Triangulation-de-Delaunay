#include <mesh.h>


#define ALPHA1 (1.0 / (4.0 * sqrt(3.0)))  // Facteur normalisation pour Q1
#define ALPHA2 (1.0 /(2.0 * sqrt(3.0)))         // Facteur normalisation pour Q2

int tri2edg[3][2] = { {1,2} , {2,0} , {0,1} };

Mesh * msh_init()
{
  Mesh *Msh = malloc(sizeof(Mesh));
  if ( ! Msh ) return NULL;
  
  Msh->Dim    = 0;
  Msh->NbrVer = 0;
  Msh->NbrTri = 0;
  Msh->NbrEfr = 0;
  Msh->NbrEdg = 0;

  Msh->Box[0] =  1.e30;   // xmin
  Msh->Box[1] = -1.e30;   // xmax
  Msh->Box[2] =  1.e30;   // ymin
  Msh->Box[3] = -1.e30;   // ymax
  
  //--- Data for the list of vertices
  Msh->Crd    = NULL;
  
  //--- Data for the list of triangles
  Msh->Tri    = NULL;
  Msh->TriVoi = NULL;
  Msh->TriRef = NULL;
  
  //--- Data for the list of boundary edges  
  Msh->Efr    = NULL;
  Msh->EfrVoi = NULL;
  Msh->EfrRef = NULL;

  //--- Data for the list of edges  
  Msh->Edg    = NULL; 
  return Msh;
}  
  
Mesh * msh_read(char *file, int readEfr)
{
  char   InpFil[1024];
  float  bufFlt[2];
  double bufDbl[2];
  int    i, bufTri[4], bufEfr[3];
  int    FilVer, ref; 
  
  int fmsh = 0;
  
  if ( ! file ) return NULL;
  
  Mesh * Msh = msh_init();
    
  //--- set file name 
  strcpy(InpFil,file);
  if ( strstr(InpFil,".mesh") ) {
    if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim)) ) {
      return NULL;
    }    
  }
  else {
    strcat(InpFil,".meshb");
    if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim)) ) {
      strcpy(InpFil,file);
      strcat(InpFil,".mesh");
      if ( !(fmsh = GmfOpenMesh(InpFil, GmfRead, &FilVer, &Msh->Dim)) ) {
        return NULL;
      }    
    } 
  }
  
  printf(" File %s opened Dimension %d Version %d \n", InpFil, Msh->Dim, FilVer);
  
  Msh->NbrVer = GmfStatKwd(fmsh, GmfVertices);
  Msh->NbrTri = GmfStatKwd(fmsh, GmfTriangles);

  Msh->NbrVerMax = Msh->NbrVer;
  Msh->NbrTriMax = Msh->NbrTri;
  
  //--- allocate arrays   
  Msh->Crd    = calloc( (Msh->NbrVer+1), sizeof(double3d) );
  Msh->Tri    = calloc( (Msh->NbrTri+1), sizeof(int3d) );  
  Msh->TriRef = calloc( (Msh->NbrTri+1), sizeof(int1d) );  
  
  
  //--- read vertices   
  GmfGotoKwd(fmsh, GmfVertices);
  if ( Msh->Dim == 2 ) {
    if ( FilVer == GmfFloat ) {		// read 32 bits float
      for (i=1; i<=Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufFlt[0], &bufFlt[1], &ref);
        Msh->Crd[i][0] = (double)bufFlt[0];
        Msh->Crd[i][1] = (double)bufFlt[1];
      }
    }
    else  {	// read 64 bits float
      for (i=1; i<=Msh->NbrVer; ++i) {
        GmfGetLin(fmsh, GmfVertices, &bufDbl[0], &bufDbl[1], &ref);
        Msh->Crd[i][0] = bufDbl[0];
        Msh->Crd[i][1] = bufDbl[1];
      }  
    }
  }
  else {
    fprintf(stderr,"  ## ERROR: 3D is not implemented\n");
    exit(1);
  }
  
  
  //--- read triangles   
  GmfGotoKwd(fmsh, GmfTriangles);
  for (i=1; i<=Msh->NbrTri; ++i) {
    GmfGetLin(fmsh, GmfTriangles, &bufTri[0], &bufTri[1], &bufTri[2], &bufTri[3]);
    Msh->Tri[i][0] = bufTri[0];
    Msh->Tri[i][1] = bufTri[1];
    Msh->Tri[i][2] = bufTri[2];
    Msh->TriRef[i] = bufTri[3];
  }
  
  
  //--- read boundary edges
  if ( readEfr == 1 ) {
    Msh->NbrEfr = GmfStatKwd(fmsh, GmfEdges);
    Msh->Efr    = calloc( (Msh->NbrEfr+1), sizeof(int2d) );  
    Msh->EfrRef = calloc( (Msh->NbrEfr+1), sizeof(int1d) );  

    GmfGotoKwd(fmsh, GmfEdges);
    for (i=1; i<=Msh->NbrEfr; ++i) {
      GmfGetLin(fmsh, GmfEdges, &bufEfr[0], &bufEfr[1], &bufEfr[2]);
      Msh->Efr[i][0] = bufEfr[0];
      Msh->Efr[i][1] = bufEfr[1];
      Msh->EfrRef[i] = bufEfr[2];
    }
  }
  
  
  GmfCloseMesh(fmsh);
  
  return Msh;
  
}

double * sol_read(char *file, int mshDim, int mshNbrSol)
{
  char   InpFil[1024];
  int    FilVer, SolTyp, NbrTyp, SolSiz, TypTab[ GmfMaxTyp ]; 
  float  bufFlt;
  double bufDbl;
  int    i, dim, nbrSol;
  
  int fsol = 0;
  
  if ( ! file ) return NULL;
  
  double * sol = NULL;
    
    
  //--- set file name 
  strcpy(InpFil, file);
  if ( strstr(InpFil,".sol") ) {
    if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
      return NULL;
    }    
  }
  else {
    strcat(InpFil,".solb");
    if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
      strcpy(InpFil,file);
      strcat(InpFil,".sol");
      if ( !(fsol = GmfOpenMesh(InpFil,GmfRead,&FilVer,&dim)) ) {
        return NULL;
      }    
    } 
  }
  
  printf(" File %s opened Dimension %d Version %d \n", InpFil, dim, FilVer);
  
  SolTyp = GmfSolAtVertices;		// read only sol at vertices
  nbrSol = GmfStatKwd(fsol, SolTyp, &NbrTyp, &SolSiz, TypTab);
	
	
  if ( nbrSol == 0 ) {
    printf("  ## WARNING: No SolAtVertices in the solution file !\n");
    return NULL;
  }
  if ( dim != mshDim ) {
    printf("  ## WARNING: WRONG DIMENSION NUMBER. IGNORED\n");
    return NULL;
  }
  if ( nbrSol != mshNbrSol ) {
    printf("  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    return NULL;
  }
  if (  NbrTyp != 1 ) {
    printf("  ## WARNING: WRONG FIELD NUMBER. IGNORED\n");
    return NULL;
  }
  if ( TypTab[0] != GmfSca ) {
    printf("  ## WARNING: WRONG FIELD TYPE. IGNORED\n");
    return NULL;
  }
	
	sol = (double *)calloc(nbrSol+1, sizeof(double));
	
	
  GmfGotoKwd(fsol, SolTyp);

  for (i=1; i<=nbrSol; ++i) {
		if ( FilVer == GmfFloat ) {
	    GmfGetLin(fsol, SolTyp, &bufFlt);
	    sol[i] = (double)bufFlt;
		}
		else {
	    GmfGetLin(fsol, SolTyp, &bufDbl);
	    sol[i] = bufDbl;
		}
  }
	
  if ( !GmfCloseMesh(fsol) ) {
    fprintf(stderr, "  ## ERROR: Cannot close solution file %s ! \n", InpFil);
    //myexit(1);
  }
	
  return sol;	
}

double area(Mesh *Msh, int iTri)
{
  int P1 = Msh->Tri[iTri][0];
  int P2 = Msh->Tri[iTri][1];
  int P3 = Msh->Tri[iTri][2];
  
  double x1 = Msh->Crd[P1][0], y1 = Msh->Crd[P1][1];
  double x2 = Msh->Crd[P2][0], y2 = Msh->Crd[P2][1];
  double x3 = Msh->Crd[P3][0], y3 = Msh->Crd[P3][1];
  
  double area = fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0;
  
  return area;
}

double inradius(Mesh *Msh, int iTri)
{
  int P1 = Msh->Tri[iTri][0];
  int P2 = Msh->Tri[iTri][1];
  int P3 = Msh->Tri[iTri][2];
  
  double x1 = Msh->Crd[P1][0], y1 = Msh->Crd[P1][1];
  double x2 = Msh->Crd[P2][0], y2 = Msh->Crd[P2][1];
  double x3 = Msh->Crd[P3][0], y3 = Msh->Crd[P3][1];
  
  double a = sqrt(pow(x2 - x3, 2) + pow(y2 - y3, 2));
  double b = sqrt(pow(x1 - x3, 2) + pow(y1 - y3, 2));
  double c = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
  
  double s = (a + b + c) / 2.0;
  double A = area(Msh, iTri);
  
  return A / s;
}

double q1_quality(Mesh *Msh, int iTri)
{
  double A = area(Msh, iTri);
  int P1 = Msh->Tri[iTri][0];
  int P2 = Msh->Tri[iTri][1];
  int P3 = Msh->Tri[iTri][2];
  
  double x1 = Msh->Crd[P1][0], y1 = Msh->Crd[P1][1];
  double x2 = Msh->Crd[P2][0], y2 = Msh->Crd[P2][1];
  double x3 = Msh->Crd[P3][0], y3 = Msh->Crd[P3][1];
  
  double a = pow(x2 - x3, 2) + pow(y2 - y3, 2);
  double b = pow(x1 - x3, 2) + pow(y1 - y3, 2);
  double c = pow(x1 - x2, 2) + pow(y1 - y2, 2);
 
  return ALPHA1 * (a+b+c) / A;
}

double q2_quality(Mesh *Msh, int iTri)
{
  double inr = inradius(Msh, iTri);
  
  int P1 = Msh->Tri[iTri][0];
  int P2 = Msh->Tri[iTri][1];
  int P3 = Msh->Tri[iTri][2];
  
  double x1 = Msh->Crd[P1][0], y1 = Msh->Crd[P1][1];
  double x2 = Msh->Crd[P2][0], y2 = Msh->Crd[P2][1];
  double x3 = Msh->Crd[P3][0], y3 = Msh->Crd[P3][1];
  
  double a = sqrt(pow(x2 - x3, 2) + pow(y2 - y3, 2));
  double b = sqrt(pow(x1 - x3, 2) + pow(y1 - y3, 2));
  double c = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
  
  double hmax = fmax(a, fmax(b, c));
  
  return ALPHA2 * (hmax / inr);
}

int msh_write(Mesh *Msh, char *file)
{
  int iVer, iTri, iEfr;
  int FilVer = 2;
  
  if ( ! Msh  ) return 0;
  if ( ! file ) return 0;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, FilVer, Msh->Dim);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  GmfSetKwd(fmsh, GmfVertices, Msh->NbrVer);
  for (iVer=1; iVer<=Msh->NbrVer; iVer++) 
    GmfSetLin(fmsh, GmfVertices, Msh->Crd[iVer][0], Msh->Crd[iVer][1], 0); 
  
  GmfSetKwd(fmsh, GmfTriangles, Msh->NbrTri);
  for (iTri=1; iTri<=Msh->NbrTri; iTri++)  
    GmfSetLin(fmsh, GmfTriangles, Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2], Msh->TriRef[iTri]);
  
  if ( Msh->NbrEfr > 0 ) {
    GmfSetKwd(fmsh, GmfEdges, Msh->NbrEfr);
    for (iEfr=1; iEfr<=Msh->NbrEfr; iEfr++)  
      GmfSetLin(fmsh, GmfEdges, Msh->Efr[iEfr][0], Msh->Efr[iEfr][1], Msh->EfrRef[iEfr]);
  }
  
  GmfCloseMesh(fmsh);
  
  return 1;
}

int msh_neighborsQ2(Mesh *Msh)
{
  int iTri, iEdg, jTri, jEdg, iVer1, iVer2, jVer1, jVer2;
  
  if ( ! Msh ) return 0;
  
  if ( Msh->TriVoi == NULL )
    Msh->TriVoi = calloc( (Msh->NbrTri+1), sizeof(int3d) );  
    for (iTri=1; iTri<=Msh->NbrTri; iTri++){
      for (iEdg=0; iEdg<3; iEdg++){
        Msh->TriVoi[iTri][iEdg]=-1;
      }

    }
  
  //--- Compute the neighbors using a quadratic-complexity algorithm 
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    if (Msh->Tri[iTri][1] != - 1){

      for (iEdg=0; iEdg<3; iEdg++) {
        iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
        iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];
        
        //--- find the Tri different from iTri that has iVer1, iVer2 as vertices 
        for (jTri=1; jTri<=Msh->NbrTri; jTri++) {
          if (( iTri == jTri ) || (Msh->Tri[jTri][1] == - 1))
            continue;
          
          for (jEdg=0; jEdg<3; jEdg++) {
            jVer1 = Msh->Tri[jTri][ tri2edg[jEdg][0] ];
            jVer2 = Msh->Tri[jTri][ tri2edg[jEdg][1] ];
            
            // TODO: compare the 4 points 
            //       set the neighbors Msh->TriVoi if both edges match
            if( ((iVer1==jVer1) && (iVer2==jVer2)) || ((iVer1==jVer2)&&(iVer2==jVer1)) ){
              Msh->TriVoi[iTri][iEdg]=jTri;
              Msh->TriVoi[jTri][jEdg]=iTri;
              break;

            }
          }
        }
        
      }
  }
}
  
  return 1;
}   

NeighborResult msh_neighbors(Mesh *Msh)
{
  int iTri, iEdg, iVer1, iVer2;
  int nbfrontiere=0;
  int nbarete=0;
  int nvoi1=0;
  int nvoi2=0;
  NeighborResult result = {0,0,0};


  FILE *file = fopen("hash_output.txt", "a"); // Ouvrir en mode ajout
    if (file == NULL) {
        perror("Erreur lors de l'ouverture du fichier");
        return result;
    }
  
  if ( ! Msh ) return result;
  
  if ( Msh->TriVoi == NULL )
    Msh->TriVoi = calloc( (Msh->NbrTri+1), sizeof(int3d) );

  
  
  //--- initialize HashTable and set the hash table 
  HashTable *hsh_tab = hash_init(2 * Msh->NbrVer, 3 * Msh->NbrTri);
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    for (iEdg=0; iEdg<3; iEdg++) {
      iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];
        hash_add(hsh_tab, iVer1, iVer2, iTri);
      }
    }
  
  for (int hashid = 1; hashid < hsh_tab->SizHead; hashid++) {
    int iobj = hsh_tab->Head[hashid];

    while (iobj != 0) { // Parcours de la liste chaînée
        fprintf(file, "num=1, hashid= %d, position du clé= %d, iobj= %d, iVer1 = %d, iVer2 = %d, iTri1 = %d, iTri2 = %d, Next= %d\n",
                hashid, hsh_tab->Head[hashid], iobj,
                hsh_tab->LstObj[iobj][0], hsh_tab->LstObj[iobj][1],
                hsh_tab->LstObj[iobj][2], hsh_tab->LstObj[iobj][3],
                hsh_tab->LstObj[iobj][4]);

        iobj = hsh_tab->LstObj[iobj][4]; // Passage à l'objet suivant dans la liste chaînée
    }
}
fclose(file);

  
  //--- Compute the neighbors using the hash table
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    for (iEdg=0; iEdg<3; iEdg++) {

      iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];

      int found = hash_find(hsh_tab, iVer1, iVer2); 
      // Récupérer les triangle déjà stocké
        int iTri1= hsh_tab->LstObj[found][2];
        int iTri2 = hsh_tab->LstObj[found][3];
        //printf("pour iVer1 = %d et iVer2 = %d , on a iTri1 = %d et iTri2 = %d \n",iVer1 , iVer2, iTri1, iTri2);
        if ((iTri2 != 0) && (iTri1 == iTri))
        {
          Msh->TriVoi[iTri][iEdg] = iTri2;
          result.nbaretevoisin +=1;
        }
        else if ((iTri1 != 0) &&(iTri2 == iTri))
        {
          Msh->TriVoi[iTri][iEdg] = iTri1;
          nvoi2+=1;

        }
        else if (((iTri2 == 0) && (iTri1 == iTri)) || ((iTri1 == 0) && (iTri2 == iTri)))
        {
          Msh->TriVoi[iTri][iEdg] = -1 ;
          result.nbfrontiere +=1 ;

        }
        }
        }
   
    printf("nvoi2 =  %d \n ", nvoi2)  ; 
    result.nbaretetotal= result.nbaretevoisin+result.nbfrontiere;
        
  return result;
}   

HashTable * hash_init(int SizHead, int NbrMaxObj)
{
	//HashTable *hsh = NULL;
	
	// to be implemented

	// allocate hash table
  HashTable *hsh = malloc(sizeof(HashTable));
	
	// initialize hash table
  hsh->NbrObj = 0;
  hsh->NbrMaxObj = NbrMaxObj;
  hsh->SizHead = SizHead;
	
	// allocate Head, LstObj
	hsh->Head = calloc(SizHead, sizeof(int));
  hsh->LstObj = calloc(NbrMaxObj, sizeof(int5d));
	
  return hsh;
}

int hash_find(HashTable *hsh, int iVer1, int iVer2)
{
  
	// to be implemented
	int hashid=hash_key(hsh, iVer1 , iVer2);
  int id=hsh->Head[hashid];
	// return the id found (in LstObj ), if 0 the object is not in the list

  while (1)
  {
    if (((hsh->LstObj[id][0] == iVer1) && (hsh->LstObj[id][1] == iVer2)) || ((hsh->LstObj[id][0] == iVer2) && (hsh->LstObj[id][1] == iVer1)))
      return id;
    if (hsh->LstObj[id][4] == 0)

      return 0;
    id = hsh->LstObj[id][4];
  }
}

int hash_key(HashTable *hsh, int iVer1 , int iVer2)
{
  return (iVer1 + iVer2) % hsh->SizHead;
   //return iVer1 >iVer2? iVer2 : iVer1;
   //return iVer1-iVer2 > 0? iVer1-iVer2 : iVer2-iVer1;
}

int hash_add(HashTable *hsh, int iVer1, int iVer2, int iTri)
{

  if (hsh->NbrObj >= hsh->NbrMaxObj - 1) return -1; // Table full
    
  int hashid = hash_key(hsh, iVer1 , iVer2);
  
  int id=hsh->Head[hashid];
    int idp=hsh->Head[hashid]; //previous on l'utilisera pour le chainage s'il existe
  if (id==0)
  {
    hsh->NbrObj +=1;
    int new_id = hsh->NbrObj;
    hsh->LstObj[new_id][0] = iVer1;
    hsh->LstObj[new_id][1] = iVer2;
    hsh->LstObj[new_id][2] = iTri;
    hsh->LstObj[new_id][3] = 0;
    hsh->LstObj[new_id][4] = 0;
    hsh->Head[hashid]=new_id;
    return 0;

  }

  if (id != 0)
  {
    // on vérifie si l'arete existe deja, dans ce cas on ajoute le iTri comme deuxieme triangle
    if (((hsh->LstObj[id][0] == iVer1) && (hsh->LstObj[id][1] == iVer2))||((hsh->LstObj[id][1] == iVer1) && (hsh->LstObj[id][0] == iVer2)))
    {
      hsh->LstObj[id][3] = iTri;
     

      return 0;
    }
    else
    {
      while (id!=0)
        {
          // clé déjà existe mais c'est new arete donc on voit la liste chainée et on garde à chaque fois la position de l'arete déjà consulté
          idp=id;  // position de l'arete previous
          id = hsh->LstObj[id][4]; // position de l'arete suivante si c'est 0, on sort du boucle while et on ajoute l'arete new
          if (((hsh->LstObj[id][0] == iVer1) && (hsh->LstObj[id][1] == iVer2))||((hsh->LstObj[id][1] == iVer1) && (hsh->LstObj[id][0] == iVer2)))
          {
            hsh->LstObj[id][3] = iTri;
            return 0;
          }
        }
        hsh->NbrObj +=1;
        int new_id = hsh->NbrObj;
        hsh->LstObj[new_id][0] = iVer1;
        hsh->LstObj[new_id][1] = iVer2;
        hsh->LstObj[new_id][2] = iTri;
        hsh->LstObj[new_id][3] = 0;
        hsh->LstObj[new_id][4] = 0;
        hsh->LstObj[idp][4] = new_id;  
    }
  }
  return 0;
}
	
int hash_suppr(HashTable *hsh, int iVer1, int iVer2, int iTri)
{
  int hashid = hash_key(hsh, iVer1 , iVer2);
  int id=hsh->Head[hashid];
  int idp=hsh->Head[hashid];
  if (id == 0){
    printf("cet objet n'existe pas");
    return 0;
  }

  while (id!=0) // ici on vérifie si la clé existe déja dans le table
  { 
    // on vérifie si l'arete existe deja, dans ce cas on ajoute le iTri comme deuxieme triangle
    if (((hsh->LstObj[id][0] == iVer1) && (hsh->LstObj[id][1] == iVer2))||((hsh->LstObj[id][1] == iVer1) && (hsh->LstObj[id][0] == iVer2)))
    {
      if (hsh->LstObj[id][3] = iTri){
        hsh->LstObj[id][3] = 0;
        printf("objet supprimé ");
        return 0;
      }
      else if (hsh->LstObj[id][2] = iTri)
      {
        hsh->LstObj[id][2] = hsh->LstObj[id][3];
        hsh->LstObj[id][3]=0;
        printf("objet supprimé ");
        break;
      }
    }
    else
    {  
      id = hsh->LstObj[id][4];
    }
  }
 if ((idp == id) && hsh->LstObj[id][2]==0)// ici cad l'objet qu'on a supprimé est la tete d'une liste chainée
 {
  hsh->Head[hashid]=hsh->LstObj[id][4];
 }

	return 0;
}

int aretefrontiere(Mesh *Msh){

  int iTri, iEdg, iVer1, iVer2;
  
  //--- initialize HashTable and set the hash table 
  HashTable *hsh_tab = hash_init(2 * Msh->NbrVer, 3 * Msh->NbrTri);
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    for (iEdg=0; iEdg<3; iEdg++) {
      iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];
        hash_add(hsh_tab, iVer1, iVer2, iTri);
      }
    }
  

  int n=0;
  for (int iObj=1 ; iObj <= hsh_tab->NbrObj ; iObj++){
    if(hsh_tab->LstObj[iObj][3] == 0) {n += 1;}
  }
  return n;
}

int collisions(Mesh *Msh){

  int iTri, iEdg, iVer1, iVer2;
  FILE *file = fopen("hash_collision.txt", "a");
  
  //--- initialize HashTable and set the hash table 
  HashTable *hsh_tab = hash_init(2 * Msh->NbrVer, 3 * Msh->NbrTri);
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    for (iEdg=0; iEdg<3; iEdg++) {
      iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];
        hash_add(hsh_tab, iVer1, iVer2, iTri);
      }
    }

    for (int hashid = 1; hashid < hsh_tab->SizHead; hashid++) {
      int iobj = hsh_tab->Head[hashid];
      int nc=1;
      while (iobj != 0) { // Parcours de la liste chaînée
        if(hsh_tab->LstObj[iobj][4] != 0) {nc += 1;}
          iobj = hsh_tab->LstObj[iobj][4]; // Passage à l'objet suivant dans la liste chaînée
      }
      if (nc > 1) {
        fprintf(file,"%d \n", nc);
      }

  }
  fclose(file);
  
  int n=0;

  for (int hashid = 1; hashid < hsh_tab->SizHead; hashid++)
   {
    int iObj = hsh_tab->Head[hashid];
    if(hsh_tab->LstObj[iObj][4] != 0) {n += 1;}
  }
  return n;
}

int msh_write2dfield_Vertices(char *file, int nfield, double *field) 
{
  int iVer;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSca;
  
  GmfSetKwd(fmsh, GmfSolAtVertices, nfield, 1, sizfld);
  
  for (iVer=1; iVer<=nfield; iVer++) 
    GmfSetLin(fmsh, GmfSolAtVertices, &field[iVer]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}

int msh_write2dfield_Triangles(char *file, int nfield, double *field) 
{
  int iTri;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSca;
  
  GmfSetKwd(fmsh, GmfSolAtTriangles, nfield, 1, sizfld);
  
  for (iTri=1; iTri<=nfield; iTri++) 
    GmfSetLin(fmsh, GmfSolAtTriangles, &field[iTri]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}

int msh_write2dmetric(char *file, int nmetric, double3d *metric) 
{  
  int iVer;
  
  int fmsh = GmfOpenMesh(file, GmfWrite, GmfDouble, 2);
  if ( fmsh <=  0 ) {
    printf("  ## ERROR: CANNOT CREATE FILE \n");
    return 0;
  }
  
  int sizfld[1];
  sizfld[0] = GmfSymMat;
  
  GmfSetKwd(fmsh, GmfSolAtVertices, nmetric, 1, sizfld);
  
  for (iVer=1; iVer<=nmetric; iVer++) 
    GmfSetLin(fmsh, GmfSolAtVertices, &metric[iVer][0], &metric[iVer][1], &metric[iVer][2]); 
  
  GmfCloseMesh(fmsh);
  
  return 1;
}




///////////////////////////////PARTiE 2/////////////////////////////////////////////////
int msh_boundingbox(Mesh *Msh)
{
  int1d iVer;
  
  //--- compute bounding box 
  for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
    // TODO: Set Msh->Box 
     if (Msh->Crd[iVer][0] <= Msh->Box[0]){
      Msh->Box[0]=Msh->Crd[iVer][0];} //xmin
    if (Msh->Crd[iVer][0] >= Msh->Box[1]){
      Msh->Box[1]=Msh->Crd[iVer][0];} //xmax
      if (Msh->Crd[iVer][1] <= Msh->Box[2]){
        Msh->Box[2]=Msh->Crd[iVer][1];} //ymin
      if (Msh->Crd[iVer][1] >= Msh->Box[3]){
        Msh->Box[3]=Msh->Crd[iVer][1]; }  //ymax
  }
  //int dx= abs(Msh->Box[1]-Msh->Box[0]);
  int dy= abs(Msh->Box[3]-Msh->Box[2]);
  return 1;
}

double areasigned(double* P1, double* P2, double* P3)
{

  double x_1 = P1[0], y_1 = P1[1];
  double x_2 = P2[0], y_2 = P2[1];
  double x_3 = P3[0], y_3 = P3[1];

  double aire = (x_1 * (y_2 - y_3) + x_2 * (y_3 - y_1) + x_3 * (y_1 - y_2)) / 2;

  return aire;
}
deuxValeurs circumcentre(double* P1, double* P2, double* P3){
  double x1 = P1[0], y1 = P1[1];
  double x2 = P2[0], y2 = P2[1];
  double x3 = P3[0], y3 = P3[1];

  deuxValeurs centre={0,0};

  double aire = areasigned(P1, P2, P3);

  centre.val1 = ((x1*x1 + y1 * y1) * y2 + (x2 * x2 + y2 * y2) * y3 + (x3 * x3 + y3 * y3) * y1 - (x3 * x3 + y3 * y3) * y2 - (x2 * x2 + y2 * y2) * y1 - (x1 * x1 + y1 * y1) * y3) / (4 * aire);
  centre.val2 = -((x1 * x1 + y1 * y1) * x2 + (x2 * x2 + y2 * y2) * x3 + (x3 * x3 + y3 * y3) * x1 - (x3 * x3 + y3 * y3) * x2 - (x2 * x2 + y2 * y2) * x1 - (x1 * x1 + y1 * y1) * x3) / (4 * aire);
  return centre;
}
troisValeurs longueur(double* P1, double* P2, double* P3)
{
  
  double x_1 = P1[0], y_1 = P1[1];
  double x_2 = P2[0], y_2 = P2[1];
  double x_3 = P3[0], y_3 = P3[1];
  troisValeurs l={0,0,0};

  l.val1 = sqrt((x_3 - x_2) * (x_3 - x_2) + (y_3 - y_2) * (y_3 - y_2));
  l.val2 = sqrt((x_3 - x_1) * (x_3 - x_1) + (y_3 - y_1) * (y_3 - y_1));
  l.val3 = sqrt((x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1));

  return l;
}
int localiser(Mesh *Msh, double* P)
{
 
  double b1, b2, b3;

  for (int iTri = 1; iTri <= Msh->NbrTri; iTri++)
  {
      if (Msh->Tri[iTri][0] > 0 && Msh->Tri[iTri][1]>0 && Msh->Tri[iTri][2]>0){
        double* P1 = Msh->Crd[Msh->Tri[iTri][0]];
        double* P2 = Msh->Crd[Msh->Tri[iTri][1]];
        double* P3 = Msh->Crd[Msh->Tri[iTri][2]];

        b1 = areasigned(P, P2, P3);
        b2 = areasigned(P1, P, P3);
        b3 = areasigned(P1, P2, P);
        // printf("%d : %f, %f, %f \n",iTri,b1,b2,b3);

        if ((b1 > 0) && (b2 >=0) && (b3 > 0))
        {
          return iTri;
        }
      }
    }
  return 0;
}

int empiler(Node **PILE, int val)
{
  Node *p = malloc(sizeof(Node));
  int success = p != NULL;

  if (success)
  {
    p->data = val;
    p->next = *PILE;
    *PILE = p;
  }

  return success;
}
int empiler_edge(NodeEdg **PILE, deuxEntiers val)
{
  NodeEdg *p = malloc(sizeof(Node));
  int success = p != NULL;

  if (success)
  {
    p->data = val;
    p->next = *PILE;
    *PILE = p;
  }

  return success;
}

int depiler(Node **PILE, int *val)
{
  int success = *PILE != NULL;

  if (success)
  {
    Node *p = *PILE;
    *PILE = (*PILE)->next;
    *val = p->data;
    free(p);
  }

  return success;
}
int depilerEdg(NodeEdg **PILE, deuxEntiers *val)
{
  int success = *PILE != NULL;

  if (success)
  {
    NodeEdg *p = *PILE;
    *PILE = (*PILE)->next;
    *val = p->data;
    free(p);
  }

  return success;
}

int InPile(Node **PILE, int val)
{

  Node *pp = *PILE;

  while (pp != NULL)
  {

    if (pp->data == val)
    {
      return 1;
    }

    pp = pp->next;
  }
  return 0;
}

int InPileEdg(NodeEdg **PILE, deuxEntiers val)
{

  NodeEdg *pp = *PILE;

  while (pp != NULL)
  {

    if ((pp->data.val1 == val.val1)&&(pp->data.val2 == val.val2))
    {
      return 1;
    }

    pp = pp->next;
  }
  return 0;
}

int VerifyInCircle(Mesh *Msh, int iTri, double *P)
{

  double* P1 = Msh->Crd[Msh->Tri[iTri][0]];
  double* P2 = Msh->Crd[Msh->Tri[iTri][1]];
  double* P3 = Msh->Crd[Msh->Tri[iTri][2]];

  troisValeurs l=longueur(P1,P2,P3);
  double aire = areasigned(P1, P2, P3);
  deuxValeurs centre=circumcentre(P1,P2,P3);

  double r_c = l.val1 * l.val2 * l.val3 / (4 * aire);

  if ((P[0] - centre.val1) * (P[0] - centre.val1) + (P[1] - centre.val2) * (P[1] - centre.val2) <= r_c * r_c)
  {
    return 1;
  }
  return 0;
}

Node *cavity(Mesh *msh, double* P)
{
  Node *pile = (NULL);
  Node *triangles=(NULL);

  int n = 0;
  int iTri = localiser(msh, P);
  
  empiler(&pile, iTri);
  //printf("triangle %d est bien dans la cavité \n" , iTri);

  for (int k = 0; k < 3; k++)
  {

    int TriVoi = msh->TriVoi[iTri][k];
    if ((TriVoi != -1) && (msh->Tri[TriVoi][0] > 0 && msh->Tri[TriVoi][1]>0 && msh->Tri[TriVoi][2]>0))
    {
      empiler(&triangles, TriVoi);
      n++;
    }
  }
  int jTri;
  while (n!=0){

    depiler(&triangles, &jTri);
    n--;
    //printf("on teste si triangle %d est dans la cavité \n", jTri);
    if (VerifyInCircle(msh, jTri, P))
    {
      empiler(&pile, jTri);
      //printf("triangle %d est bien dans la cavité \n" , jTri);

      for (int k = 0; k < 3; k++)
      {
        int VoiK = msh->TriVoi[jTri][k];
        if ((VoiK != -1) && !(InPile(&triangles, VoiK)) && !(InPile(&pile, VoiK))&& (msh->Tri[VoiK][1] != -1)) 
        {
          empiler(&triangles, VoiK);
          //printf("on ajoute  %d est dans triangles à tester \n", VoiK);
          n++;
        }
      }
    }
  }
  
  return pile;
}

void print_pile(Node *pile)
{
  Node *pp = pile;
  while (pp != NULL)
  {  
    printf("%d  \n ", pp->data);
    pp = pp->next;
  }
}
void print_pile_Edge(NodeEdg *pile)
{
  NodeEdg *pp = pile;
  while (pp != NULL)
  {  
    printf("iedg 1 : %d \n ", pp->data.val1);
    printf("iedg 2 : %d \n ", pp->data.val2);
    pp = pp->next;
  }
}

NodeEdg *frontiere(Mesh *msh, Node* cavité){
  Node *cavité_temp = cavité;
  if (cavité==NULL){
    printf("erreur, cavité vide");
    return NULL ;
  }
  NodeEdg *pileEdg = (NULL);

  while (cavité_temp != NULL) {
      int iTri = cavité_temp->data;

      for (int iEdg = 0; iEdg < 3; iEdg++) {
          int voisin = msh->TriVoi[iTri][iEdg];
          
              // Si l'arête est en bordure de la cavité
              if (voisin == -1 || InPile(&cavité, voisin)==0 ) {
                deuxEntiers edg;
                  edg.val1=msh->Tri[iTri][ tri2edg[iEdg][0] ];  // Premier sommet de l'arête
                  edg.val2=msh->Tri[iTri][ tri2edg[iEdg][1] ];  // Deuxième sommet de l'arête
                  empiler_edge(&pileEdg,edg);
                  //printf("pour le triangle %d on note %d arete de %d et %d \n",iTri,iEdg,edg.val1,edg.val2);   
              }
          
    }

      cavité_temp = cavité_temp->next;
  }
  return pileEdg;
}
int efface_contenu(Mesh *Msh, double *P){
  int iTri = localiser(Msh, P);
  if (iTri ==0){
    printf("le point (%f ,%f) n'existe pas dans le maillage \n",P[0],P[1]);
    return 0;
  }

  //on crée la cavité et sa frontiere
  Node *cavité = cavity(Msh, P);
  Node *cavité_temp=cavité;
  Node *cavité_ref=cavité;
  NodeEdg *CavityEdg =frontiere(Msh,cavité);
  NodeEdg *CavityEdg_temp =CavityEdg;



  while (cavité != NULL) {
    int jTri = cavité->data;  // Récupération de l'indice du triangle
        Msh->Tri[jTri][0] = - Msh->Tri[jTri][0];
    // Passer au prochain élément de la cavité
    cavité = cavité->next;
}
}

int delaunay(Mesh *Msh, double *P)
{ 
  //localisation du point déjà dans la fonction cavity
  int iTri = localiser(Msh, P);
  if (iTri ==0){
    printf("le point (%f ,%f) n'existe pas dans le maillage \n",P[0],P[1]);
    return 0;
  }

  //on crée la cavité et sa frontiere
  Node *cavité = cavity(Msh, P);
  Node *cavité_temp=cavité;
  NodeEdg *CavityEdg =frontiere(Msh,cavité);
  NodeEdg *CavityEdg_temp =CavityEdg;
  int newtriangle=0;

  //printf("NberVer= %d \n" , Msh->NbrVer);
  Msh->NbrVer++;

  Msh->Crd = realloc(Msh->Crd, (Msh->NbrVer + 3) * sizeof(double2d));
  Msh->NbrVerMax *= 2;
  
  //printf("new NberVer= %d \n" , Msh->NbrVer);
  Msh->Crd[Msh->NbrVer][0] = P[0];
  Msh->Crd[Msh->NbrVer][1] = P[1];
  //printf("ajout de P \n ");

  while (cavité != NULL) {
    int jTri = cavité->data;  // Récupération de l'indice du triangle
        Msh->Tri[jTri][1] = - 1; //mettre le triangle inactif en mettant -1 sur l'un 
    cavité = cavité->next;
}

  while (CavityEdg_temp != (NULL)){
    Msh->NbrTri++;

    Msh->Tri = realloc(Msh->Tri, (Msh->NbrTri + 1) * sizeof(int3d));
    Msh->TriVoi = realloc(Msh->TriVoi, (Msh->NbrTri + 1) * sizeof(int3d));
    Msh->TriRef = realloc(Msh->TriRef, (Msh->NbrTri + 1) * sizeof(int1d));
          
    
    int iTri=Msh->NbrTri;
    deuxEntiers edg = CavityEdg_temp->data;
    Msh->Tri[iTri][0]=Msh->NbrVer;
    Msh->Tri[iTri][1]=edg.val1 > 0? edg.val1 : -edg.val1;
    Msh->Tri[iTri][2]=edg.val2>0 ? edg.val2 : -edg.val2;
    newtriangle++;




    //maj les voisins
    // Initialisation des voisins à -1
    Msh->TriVoi[iTri][0] = -1;
    Msh->TriVoi[iTri][1] = -1;
    Msh->TriVoi[iTri][2] = -1;

    Msh->TriRef[iTri] = 1;

    for (int iEdg=0; iEdg<3; iEdg++) {
      int iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      int iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];
      
      //--- find the Tri different from iTri that has iVer1, iVer2 as vertices 
      for (int jTri=1; jTri< Msh->NbrTri; jTri++) {
        if (Msh->Tri[jTri][1] != - 1){

        for (int jEdg=0; jEdg<3; jEdg++) {
          int jVer1 = Msh->Tri[jTri][ tri2edg[jEdg][0] ];
          int jVer2 = Msh->Tri[jTri][ tri2edg[jEdg][1] ];
          
          // TODO: compare the 4 points 
          //       set the neighbors Msh->TriVoi if both edges match
          if( ((iVer1==jVer1) && (iVer2==jVer2)) || ((iVer1==jVer2)&&(iVer2==jVer1)) ){
            Msh->TriVoi[iTri][iEdg]=jTri;
            Msh->TriVoi[jTri][jEdg]=iTri;
            //printf("%d est le voisin de %d par %d - %d \n", iTri,jTri,iVer1,iVer2);
            break;

          }
        }
      }
      }
      
    }

    CavityEdg_temp=CavityEdg_temp->next;
  }
  //printf("%d triangles ajoutés \n", newtriangle);
  return newtriangle;
}

int msh_split_triangle(Mesh *Msh, int iTri, double *P) {
  if (iTri < 0 || iTri >= Msh->NbrTri) {
      printf("Erreur: Indice de triangle invalide.\n");
      return -1;
  }
  
  // Récupération des sommets du triangle à diviser
  int iVer1 = Msh->Tri[iTri][0];
  int iVer2 = Msh->Tri[iTri][1];
  int iVer3 = Msh->Tri[iTri][2];
  
  // Ajout du nouveau sommet aux coordonnées
  int NewIVer = ++Msh->NbrVer;
  Msh->Crd = realloc(Msh->Crd, (Msh->NbrVer + 1) * sizeof(double2d));
  Msh->Crd[NewIVer][0] = P[0];
  Msh->Crd[NewIVer][1] = P[1];
 
  
  // Création des trois nouveaux triangles
  int newTriIdx1 = iTri; // Remplace l'ancien triangle
  int newTriIdx2 = Msh->NbrTri+1;
  int newTriIdx3 = Msh->NbrTri+2;
  Msh->Tri = realloc(Msh->Tri, (Msh->NbrTri + 3) * sizeof(int3d));
  Msh->TriVoi = realloc(Msh->TriVoi, (Msh->NbrTri + 3) * sizeof(int3d));
  Msh->TriRef = realloc(Msh->TriRef, (Msh->NbrTri + 3) * sizeof(int1d));
  
  // Mise à jour des nouveaux triangles
  Msh->Tri[newTriIdx1][0] = iVer1; Msh->Tri[newTriIdx1][1] = iVer2; Msh->Tri[newTriIdx1][2] = NewIVer;
  Msh->Tri[newTriIdx2][0] = iVer2; Msh->Tri[newTriIdx2][1] = iVer3; Msh->Tri[newTriIdx2][2] = NewIVer;
  Msh->Tri[newTriIdx3][0] = iVer3; Msh->Tri[newTriIdx3][1] = iVer1; Msh->Tri[newTriIdx3][2] = NewIVer;
  
  // Mise à jour des références
  Msh->TriRef[newTriIdx1] = Msh->TriRef[iTri];
  Msh->TriRef[newTriIdx2] = Msh->TriRef[iTri];
  Msh->TriRef[newTriIdx3] = Msh->TriRef[iTri];
  
  // Mise à jour des voisins
  Msh->TriVoi[newTriIdx1][0] = newTriIdx2; Msh->TriVoi[newTriIdx1][1] = newTriIdx3; Msh->TriVoi[newTriIdx1][2] = Msh->TriVoi[iTri][0];
  Msh->TriVoi[newTriIdx2][0] = newTriIdx3; Msh->TriVoi[newTriIdx2][1] = newTriIdx1; Msh->TriVoi[newTriIdx2][2] = Msh->TriVoi[iTri][1];
  Msh->TriVoi[newTriIdx3][0] = newTriIdx1; Msh->TriVoi[newTriIdx3][1] = newTriIdx2; Msh->TriVoi[newTriIdx3][2] = Msh->TriVoi[iTri][2];
  Msh->NbrTri=Msh->NbrTri+2;
  return newTriIdx1;
}
void detect_contours_sobel(Mesh* msh, double* sol, double seuil, double2d *ListeP, double* Sol_Com, int* n) {
  msh_boundingbox(msh);
  int Nx = msh->Box[1];
  int Ny = msh->Box[3]; 

  for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
          // Calcul de l'indice iVer du point (i, j)
          int idx = i + (Ny-j) * Nx;

          // Appliquer les noyaux Sobel
          double gx =
              -sol[(i - 1) + (Ny-(j - 1)) * Nx] - 2 * sol[(i - 1) + (Ny-j)* Nx] - sol[(i - 1) + (Ny-(j + 1)) * Nx]
              + sol[(i + 1) +(Ny-(j - 1)) * Nx] + 2 * sol[(i + 1) + (Ny-j)* Nx] + sol[(i + 1) + (Ny-(j + 1)) * Nx];

          double gy =
              -sol[(i - 1) + (Ny-(j - 1)) * Nx] - 2 * sol[i +(Ny-(j - 1))* Nx] - sol[(i + 1) + (Ny-(j - 1)) * Nx]
              + sol[(i - 1) + (Ny- (j + 1)) * Nx] + 2 * sol[i + (Ny-(j + 1)) * Nx] + sol[(i + 1) + (Ny-(j + 1)) * Nx];

          double grad = sqrt(gx * gx + gy * gy);

          if (grad > seuil) {
              int iVer = idx;
              ListeP[*n][0] = msh->Crd[iVer][0];
              ListeP[*n][1] = msh->Crd[iVer][1];
              Sol_Com[*n + 4] = sol[iVer];
              (*n)++;
          }
      }
  }
}

double PSNR(Mesh *msh_compressed, Mesh *msh_origine, double* sol_compressed, double* sol_originale) {
  msh_boundingbox(msh_origine);
  int Nx = msh_origine->Box[1];  // nombre de points en x
  int Ny = msh_origine->Box[3];  // nombre de points en y
  int jmax = Ny;                 

  double qc = 0.0;
  int count = 0;
  double d = 255.0;

  for (int j = 1; j <= Ny; j++) {
      for (int i = 1; i <= Nx; i++) {
          int iVer = i + (jmax - j) * Nx;

          double P[2] = { msh_origine->Crd[iVer][0], msh_origine->Crd[iVer][1] };
          int iTri = localiser(msh_compressed, P);

          if (iTri >= 0) {
              int* tri = msh_compressed->Tri[iTri];
              double* P1 = msh_compressed->Crd[tri[0]];
              double* P2 = msh_compressed->Crd[tri[1]];
              double* P3 = msh_compressed->Crd[tri[2]];

              double beta1 = areasigned(P, P2, P3) / areasigned(P1, P2, P3);
              double beta2 = areasigned(P1, P, P3) / areasigned(P1, P2, P3);
              double beta3 = areasigned(P1, P2, P) / areasigned(P1, P2, P3);

              double uc = beta1 * sol_compressed[tri[0]] +
                          beta2 * sol_compressed[tri[1]] +
                          beta3 * sol_compressed[tri[2]];

              double u = sol_originale[iVer];
              qc += (u - uc) * (u - uc);
              count++;
          }
      }
  }

  if (count == 0) {
      fprintf(stderr, "Erreur : aucun point localisé dans le maillage compressé.\n");
      return -1;
  }

  qc /= count;
  double psnr = 10.0 * log10((d * d) / qc);
  return psnr;
}


