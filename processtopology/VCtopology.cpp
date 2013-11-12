/***************************************************************************
                  VCtopology.h  -  Virtual cartesian topology
                  A virtual topology is a mechanism for naming the processes
                  in a communicator in a way that fits the communication
                  pattern better. Since our processes will communicate mainly
                  with the nearest neighbours after the fashion of a two-dimensional
                  grid, we create a virtual topology to reflect this fact
                             -------------------

 ***************************************************************************/


#include "mpi.h"

#include <iostream>

#include "VCtopology.h"

using std::cout;
using std::endl;
using std::cerr;
/**
*
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta.
* @version 2.0
*
*/

/** constructor: initialize a virtual cartesian topology */
VCtopology::VCtopology(){
    /*************************************************************/
    /***  HERE YOU SET THE TOPOLOGY    XLEN*YLEN = n_procs     ***/
    /*************************************************************/
    /****** PERIODICITY OF PROCESS TOPOLOGY   ***/
    /** IF ONE OF THIS PERIODICITY IS TRUE PERIODICITY IS AUTOMATICALLY IMPOSED ON BC !!!!*/
    /** Periodicity in Z must always be false since Z is the grid index **/
    /**************************************************************/

    PERIODICX = false;
    PERIODICY = false;
    XDIR = 0;
    YDIR = 1;
    RIGHT = 1;
    LEFT = -1;
    reorder = 1;
    divisions = new int[3];
    periodic_divisions = new int[2];
    periods = new int[3];
    periods[0] = false;
    periods[1] = false;
    periods[2] = false;
    coordinates = new int[2];
    global_coordinates = new int[3];
    cVERBOSE = false;
    twoDperiods = new int[2];

}
/** destructor */
VCtopology::~VCtopology(){
  delete[] periods;
  delete[] divisions;
  delete[] coordinates;
  delete[] global_coordinates;
  delete[] periodic_divisions;
  delete[] twoDperiods;
}
/** Within CART_COMM, processes find about their new rank numbers, their cartesian coordinates,
    and their neighbors  */
void VCtopology::setup_vctopology(MPI_Comm old_comm){
      int remain_dims[3];
      PROCDIM =2;
      remain_dims[0]= 0;
      remain_dims[1]= 1;      
      remain_dims[2]= 1;      

    MPI_Cart_create(old_comm, 3, divisions, periods, reorder, &CART_COMM_TOTAL);
    MPI_Cart_sub(CART_COMM_TOTAL,remain_dims, &CART_COMM_test);
    // Trying to make the coarsest level communicator periodic
    MPI_Comm_rank   (CART_COMM_TOTAL, &cartesian_rank);
    MPI_Cart_coords(CART_COMM_TOTAL, cartesian_rank, 3, global_coordinates);
    if (global_coordinates[0] == 0){  //If level == 0
        PERIODICX = twoDperiods[0];
        PERIODICY = twoDperiods[1];
        if (cartesian_rank == 0){
            if (PERIODICX)
                cout << "Topology is periodic in the X direction"<<endl;
            if (PERIODICY)
                cout << "Topology is periodic in the Y direction"<<endl;
        }
        periodic_divisions[0] = divisions[1];
        periodic_divisions[1] = divisions[2];
        MPI_Cart_create(CART_COMM_test, 2, periodic_divisions, twoDperiods, reorder, &CART_COMM);
    } else { 
    CART_COMM = CART_COMM_test;
    }

    if (CART_COMM != MPI_COMM_NULL){
      int *tempCoor = new int[2];
      int bx, by;
      MPI_Comm_rank   (CART_COMM, &cartesian_rank);
      MPI_Cart_coords (CART_COMM, cartesian_rank, 2, coordinates);


      MPI_Cart_shift  (CART_COMM, XDIR, RIGHT, &xleft_neighbor, &xright_neighbor);
      MPI_Cart_shift  (CART_COMM, YDIR, RIGHT, &yleft_neighbor, &yright_neighbor);
      //cout << cartesian_rank<<" "<<xleft_neighbor << " " << xright_neighbor<<" " <<yleft_neighbor<<" "<<yright_neighbor<<endl;
      // know the 12 processors in diagonal that share edges
      /** cartesian rank of XLEFT(-) YLEFT(-) */
      if ( (xleft_neighbor!= MPI_PROC_NULL) && (yleft_neighbor!= MPI_PROC_NULL) ){
          bx =   ((coordinates[0])==(0))?(XLEN):(0);
          by =   ((coordinates[1])==(0))?(YLEN):(0);

          tempCoor[0] = (coordinates[0]-1) + bx;
          tempCoor[1] = (coordinates[1]-1) + by;
          MPI_Cart_rank(CART_COMM,tempCoor,&XleftYleft_neighbor);

      } else {
           XleftYleft_neighbor = MPI_PROC_NULL;

      }
      /** cartesian rank of XLEFT(-) YRIGHT(+)*/
      if( (xleft_neighbor!= MPI_PROC_NULL) && (yright_neighbor!= MPI_PROC_NULL) ){
        bx = ((coordinates[0])==(0))?(XLEN):(0);
        tempCoor[0] = (coordinates[0]-1) + bx;
        tempCoor[1] = (coordinates[1]+1)%(YLEN);
        MPI_Cart_rank(CART_COMM,tempCoor,&XleftYright_neighbor);

      } else
      {
         XleftYright_neighbor =  MPI_PROC_NULL;
      }
      /** cartesian rank of XRIGHT(+) YLEFT(-) SAME Z neighbor */
      if((xright_neighbor!= MPI_PROC_NULL) && (yleft_neighbor!= MPI_PROC_NULL)){
        by =  ((coordinates[1])==(0))?(YLEN):(0);
        tempCoor[0] = (coordinates[0]+1)%XLEN;
        tempCoor[1] = coordinates[1]-1 + by;
        MPI_Cart_rank(CART_COMM,tempCoor,&XrightYleft_neighbor);
      }  else
      {
         XrightYleft_neighbor =  MPI_PROC_NULL;
      }
      /** cartesian rank of XRIGHT(+) YRIGHT(+) SAME Z neighbor */
      if ((xright_neighbor!= MPI_PROC_NULL) && (yright_neighbor!= MPI_PROC_NULL)){
        tempCoor[0] = (coordinates[0]+1)%XLEN;
        tempCoor[1] = (coordinates[1]+1)%YLEN;
        MPI_Cart_rank(CART_COMM,tempCoor,&XrightYright_neighbor);
      } else
      {
        XrightYright_neighbor =  MPI_PROC_NULL;
      }


    } else
    {
      cout << "A process is trown away from the new topology. VCtopology.h" << endl;
    }

    // set up a communicator per direction for refined level boundary processors
    int* ranks_leftB= new int[YLEN];
    int* ranks_rightB= new int[YLEN];
    int* ranks_bottomB= new int[XLEN];
    int* ranks_topB= new int[XLEN];
    int* ranks_AllB= new int[2*(XLEN+YLEN-2)];

    // rank on the boundary communicators
    COMM_B_LEFT=MPI_COMM_NULL, COMM_B_RIGHT=MPI_COMM_NULL, COMM_B_BOTTOM=MPI_COMM_NULL, COMM_B_TOP=MPI_COMM_NULL, COMM_B_ALL=MPI_COMM_NULL;

    int rankOnLevel= getCartesian_rank();

    MPI_Group orig_group, new_groupL, new_groupR, new_groupB, new_groupT, new_groupAllB;

    if (global_coordinates[0]==1)   //level 1
      {
	int i_left=0, i_right=0, i_bottom=0, i_top=0, i_AllB=0;
	int numtasks;
	
	MPI_Comm_group(CART_COMM, &orig_group);
	MPI_Comm_size(CART_COMM, &numtasks);
	for (int i=0; i<numtasks; i++) // create the list of members for each group
	  {
	    if(i<=YLEN-1)
	      {
		ranks_leftB[i_left]=i; i_left++;
		if (i!=0 and i!=YLEN-1)
		  {ranks_AllB[i_AllB]=i; i_AllB++;}
	      }
	    if(i>=numtasks-YLEN)
	      {
		ranks_rightB[i_right]=i; i_right++;
		if (i!=numtasks-YLEN and i!=XLEN*YLEN-1)
		  {ranks_AllB[i_AllB]=i; i_AllB++;}
	      }
	    if(i%YLEN==0)
	      {ranks_bottomB[i_bottom]=i; i_bottom++;
		ranks_AllB[i_AllB]=i; i_AllB++;}
	    if(i%YLEN==YLEN-1)
	      {ranks_topB[i_top]=i; i_top++;
		ranks_AllB[i_AllB]=i; i_AllB++;}
	  }


	if (i_left!=YLEN or i_right!=YLEN or i_bottom!=XLEN or i_top!=XLEN or i_AllB!=2*(XLEN+YLEN-2))
	  {
	    cerr <<"ERROR IN THE CREATION OF THE BOUNDARY COMMUNICATORS!"<<endl;
	    cout <<"ERROR IN THE CREATION OF THE BOUNDARY COMMUNICATORS!"<<endl;
	    return;
	  }

	int returnL, returnR, returnB, returnT;

	MPI_Group_incl(orig_group, YLEN, ranks_leftB, &new_groupL);
	MPI_Group_incl(orig_group, YLEN, ranks_rightB, &new_groupR);
	MPI_Group_incl(orig_group, XLEN, ranks_bottomB, &new_groupB);
	MPI_Group_incl(orig_group, XLEN, ranks_topB, &new_groupT);

	MPI_Group_incl(orig_group, 2*(XLEN+YLEN-2), ranks_AllB, &new_groupAllB);
	  
	returnL=MPI_Comm_create(CART_COMM, new_groupL, &COMM_B_LEFT); 
	returnR=MPI_Comm_create(CART_COMM, new_groupR, &COMM_B_RIGHT);
	returnB=MPI_Comm_create(CART_COMM, new_groupB, &COMM_B_BOTTOM);
	returnT=MPI_Comm_create(CART_COMM, new_groupT, &COMM_B_TOP);

	MPI_Comm_create(CART_COMM, new_groupAllB, &COMM_B_ALL);

	neighborLeft_COMM_B_LEFT=MPI_PROC_NULL, neighborRight_COMM_B_LEFT=MPI_PROC_NULL;
	rankOnBLeft=-1;
	//neighbors and rank
	if (returnL== MPI_SUCCESS and COMM_B_LEFT!= MPI_COMM_NULL)   
	  {
	    MPI_Comm_rank(COMM_B_LEFT, &rankOnBLeft);
	    if (rankOnBLeft>0) {neighborLeft_COMM_B_LEFT=rankOnBLeft-1;}
	    else {neighborLeft_COMM_B_LEFT=MPI_PROC_NULL;}
	    if (rankOnBLeft<YLEN-1){neighborRight_COMM_B_LEFT=rankOnBLeft+1;}
	    else {neighborRight_COMM_B_LEFT=MPI_PROC_NULL;}
	  }

	neighborLeft_COMM_B_RIGHT=MPI_PROC_NULL, neighborRight_COMM_B_RIGHT=MPI_PROC_NULL;
        rankOnBRight=-1;
        //neighbors and rank                                                           
        if (returnR== MPI_SUCCESS and COMM_B_RIGHT!= MPI_COMM_NULL)
          {
            MPI_Comm_rank(COMM_B_RIGHT, &rankOnBRight);
	    if (rankOnBRight>0) {neighborLeft_COMM_B_RIGHT=rankOnBRight-1;}
            else {neighborLeft_COMM_B_RIGHT=MPI_PROC_NULL;}
            if (rankOnBRight<YLEN-1){neighborRight_COMM_B_RIGHT=rankOnBRight+1;}
	    else {neighborRight_COMM_B_RIGHT=MPI_PROC_NULL;}
          }

	neighborLeft_COMM_B_BOTTOM=MPI_PROC_NULL, neighborRight_COMM_B_BOTTOM=MPI_PROC_NULL;
        rankOnBBottom=-1;
        //neighbors and rank                                   
        if (returnB== MPI_SUCCESS and COMM_B_BOTTOM!= MPI_COMM_NULL)
          {
            MPI_Comm_rank(COMM_B_BOTTOM, &rankOnBBottom);
	    if (rankOnBBottom>0) {neighborLeft_COMM_B_BOTTOM=rankOnBBottom-1;}
            else {neighborLeft_COMM_B_BOTTOM=MPI_PROC_NULL;}
            if (rankOnBBottom<XLEN-1){neighborRight_COMM_B_BOTTOM=rankOnBBottom+1;}
	    else {neighborRight_COMM_B_BOTTOM=MPI_PROC_NULL;}
          }
	
	neighborLeft_COMM_B_TOP=MPI_PROC_NULL, neighborRight_COMM_B_TOP=MPI_PROC_NULL;
        rankOnBTop=-1;
        //neighbors and rank                                                            
        if (returnT== MPI_SUCCESS and COMM_B_TOP!= MPI_COMM_NULL)
          {
            MPI_Comm_rank(COMM_B_TOP, &rankOnBTop);
	    if (rankOnBTop>0) {neighborLeft_COMM_B_TOP=rankOnBTop-1;}
            else {neighborLeft_COMM_B_TOP=MPI_PROC_NULL;}
            if (rankOnBTop<XLEN-1){neighborRight_COMM_B_TOP=rankOnBTop+1;}
            else {neighborRight_COMM_B_TOP=MPI_PROC_NULL;}
          }
	
	/*//DEBUG
	MPI_Barrier(CART_COMM);                                                                                     
        if (rankOnLevel==0)                                                           
          {                                                                                   
	    cout <<"ranks_AllB is :\n";                                                        
	    for (int i=0; i<i_AllB; i++) cout<<ranks_AllB[i] <<endl;  
	  }
	MPI_Barrier(CART_COMM);
	MPI_Barrier(CART_COMM);
        if (rankOnLevel==0)
          {
            cout <<"ranks_leftB is :\n";
            for (int i=0; i<i_left; i++) cout<<ranks_leftB[i] <<endl;
            cout <<"ranks_rightB is :\n";
            for(int i=0; i<i_right; i++) cout<<ranks_rightB[i] <<endl;
            cout <<"ranks_bottomB is :\n";
            for(int i=0; i<i_bottom; i++) cout<<ranks_bottomB[i] <<endl;
            cout <<"ranks_topB is :\n";
            for(int i=0; i<i_top; i++) cout<<ranks_topB[i] <<endl;
	    
          }
        MPI_Barrier(CART_COMM);
     
	if (returnL== MPI_SUCCESS and COMM_B_LEFT!= MPI_COMM_NULL)
	  {
	    //rankOnBLeft=getRank_COMM_B_LEFT();
	    cout << "My level rank is "<<rankOnLevel <<" and I am part of COMM_B_LEFT\n";}
	if (returnR== MPI_SUCCESS and COMM_B_RIGHT!= MPI_COMM_NULL)
	  {
	    //rankOnBRight=getRank_COMM_B_RIGHT();
	    cout << "My level rank is "<<rankOnLevel <<" and I am part of COMM_B_RIGHT\n";}
	if (returnB== MPI_SUCCESS and COMM_B_BOTTOM!= MPI_COMM_NULL)
	  {
	    //rankOnBBottom=getRank_COMM_B_BOTTOM();
	    cout << "My level rank is "<<rankOnLevel <<" and I am part of COMM_B_BOTTOM\n";}
	if (returnT== MPI_SUCCESS and COMM_B_TOP!= MPI_COMM_NULL)
	  {
	    //rankOnBTop=getRank_COMM_B_TOP();
	    cout << "My level rank is "<<rankOnLevel <<" and I am part of COMM_B_TOP\n";}

	int testlocalL=1, testglobalL=0, testlocalR=1, testglobalR=0, testlocalB=1, testglobalB=0, testlocalT=1, testglobalT=0;
	if (returnL== MPI_SUCCESS and COMM_B_LEFT!= MPI_COMM_NULL)
	  MPI_Reduce(&testlocalL, &testglobalL, 1, MPI_INT, MPI_SUM, 0, COMM_B_LEFT);
	if (returnR== MPI_SUCCESS and COMM_B_RIGHT!= MPI_COMM_NULL)
	  MPI_Reduce(&testlocalR, &testglobalR, 1, MPI_INT, MPI_SUM, 0, COMM_B_RIGHT);
	if (returnB== MPI_SUCCESS and COMM_B_BOTTOM!= MPI_COMM_NULL)
	  MPI_Reduce(&testlocalB, &testglobalB, 1, MPI_INT, MPI_SUM, 0, COMM_B_BOTTOM);
	if (returnT== MPI_SUCCESS and COMM_B_TOP!= MPI_COMM_NULL)
	  MPI_Reduce(&testlocalT, &testglobalT, 1, MPI_INT, MPI_SUM, 0, COMM_B_TOP);

	if (rankOnLevel==ranks_leftB[0])
	  {cout << testglobalL <<" processes are part of the communicator COMM_B_LEFT, they should be " << YLEN <<endl; }
	if (rankOnLevel==ranks_rightB[0])
          {cout << testglobalR <<" processes are part of the communicator COMM_B_RIGHT, they should be " << YLEN <<endl; }
	if (rankOnLevel==ranks_bottomB[0])
          {cout << testglobalB <<" processes are part of the communicator COMM_B_BOTTOM, they should be " << XLEN <<endl; }
	if (rankOnLevel==ranks_topB[0])
          {cout << testglobalT <<" processes are part of the communicator COMM_B_TOP, they should be " << XLEN <<endl; }  
	  // End DEBUG*/
	
      } // end level 1

    delete[] ranks_leftB;
    delete[] ranks_rightB;
    delete[] ranks_bottomB;
    delete[] ranks_topB;
    delete[] ranks_AllB;
}

/** print topology info */
void VCtopology::Print(){
  cout << endl;
  cout << "Virtual Cartesian Processors Topology" << endl;
  cout << "-------------------------------------"<< endl;
  cout << "Processors Topology dimension = " << PROCDIM << endl;
  cout << "Processors grid: " << XLEN << "x" << YLEN << endl;
  cout << "Periodicity X: " << PERIODICX << endl;
  cout << "Periodicity Y: " << PERIODICY << endl;
  cout << "Number of levels : " << ngrids << endl;
  cout << endl;
}
/** print cartesian rank of neighbors and coordinate of process */
void VCtopology::PrintMapping(){
  cout << endl;
  cout << "Mapping of process " << cartesian_rank << endl;
  cout << "----------------------" << endl;
  cout << "Coordinates: X = " << coordinates[0] << "; Y = " << coordinates[1] << endl;
  cout << "Neighbors: xLeft = " << xleft_neighbor << "; xRight = " << xright_neighbor << "; yLeft = " << yleft_neighbor << "; yRight = "<< yright_neighbor << endl;
  cout << "Neighbors: xRightYright = " << XrightYright_neighbor << "; xRightYleft = " << XrightYleft_neighbor << "; xLeftYright = " << XleftYright_neighbor << "; xLeftYleft = "<< XleftYleft_neighbor << endl;

  cout << endl;
}
/** get and set XLEN */
int VCtopology::getXLEN(){
  return(XLEN);
}
void VCtopology::setXLEN(int XLENinput){
  XLEN = XLENinput;
}
/** get and set YLEN */
int VCtopology::getYLEN(){
  return(YLEN);
}
void VCtopology::setYLEN(int YLENinput){
  YLEN = YLENinput;
}
/** get CART_COMM and CART_COMM_TOTAL**/
MPI_Comm VCtopology::getCART_COMM(){
   return(CART_COMM);
}
MPI_Comm VCtopology::getCART_COMM_TOTAL(){
   return(CART_COMM_TOTAL);
}
/** get boundary communicators **/
MPI_Comm VCtopology::getCOMM_B_LEFT(){
  return(COMM_B_LEFT);
}
MPI_Comm VCtopology::getCOMM_B_RIGHT(){
  return(COMM_B_RIGHT);
}
MPI_Comm VCtopology::getCOMM_B_BOTTOM(){
  return(COMM_B_BOTTOM);
}
MPI_Comm VCtopology::getCOMM_B_TOP(){
  return(COMM_B_TOP);
}
MPI_Comm VCtopology::getCOMM_B_ALL(){
  return(COMM_B_ALL);
}
/** get and set ngrids */
int VCtopology::getNgrids(){
  return(ngrids);
}
void VCtopology::setNgrids(int ngridinput){
  ngrids = ngridinput;
}
/** set divisions **/
void VCtopology::setDivisions(){
    divisions[0] = ngrids;
    divisions[1] = XLEN;
    divisions[2] = YLEN;
}

/** get and set nprocs */
int VCtopology::getNprocs(){
  return(nprocs);
}
void VCtopology::setNprocs(){
  //cout << "XLEN*YLEN= " << XLEN*YLEN << endl;
  nprocs = XLEN*YLEN;
}
/** get periodicity on boundaries - DIRECTION X*/
bool VCtopology::getPERIODICX(){
  return(PERIODICX);
}
/** get periodicity on boundaries - DIRECTION Y*/
bool VCtopology::getPERIODICY(){
  return(PERIODICY);
}
/** set periodicity **/
void VCtopology::setPeriodicity(int xleft,int xright,int yleft,int yright){
  if(xleft == 2|| xright ==2){
      twoDperiods[0] = true;
  } else {
     twoDperiods[0] = false;
  }
  if(yleft == 2|| yright ==2){
      twoDperiods[1] = true;
  } else {
      twoDperiods[1] = false;
  }
}
/** get the cartesian rank of the process */
int VCtopology::getCartesian_rank(){
  return(cartesian_rank);
}
/** get the rank of the process in MPI_COMM_TOTAL */                                                      int VCtopology::getCartesian_rank_COMMTOTAL(){
  int r;
  MPI_Comm_rank(CART_COMM_TOTAL, &r);
  return r;
}
int VCtopology::getRank_MPI_COMM_WORLD(){
  int r;
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  return r;
}
/** get the rank of the process on the boundary communicators */
int VCtopology::getRank_COMM_B_LEFT(){
  /*int r;
  MPI_Comm_rank(COMM_B_LEFT, &r);
  return r;*/
  return rankOnBLeft;
}
int VCtopology::getRank_COMM_B_RIGHT(){
  /*int r;
  MPI_Comm_rank(COMM_B_RIGHT, &r);
  return r;*/
  return rankOnBRight;
}
int VCtopology::getRank_COMM_B_BOTTOM(){
  /*int r;
  MPI_Comm_rank(COMM_B_BOTTOM, &r);
  return r;*/
  return rankOnBBottom;
}
int VCtopology::getRank_COMM_B_TOP(){
  /*int r;
  MPI_Comm_rank(COMM_B_TOP, &r);
  return r;*/
  return rankOnBTop;
}
/** get the cartesian rank of XLEFT neighbor */
int VCtopology::getXleft_neighbor(){
  return(xleft_neighbor);
}
/** get the cartesian rank of XRIGHT neighbor */
int VCtopology::getXright_neighbor(){
  return(xright_neighbor);
}
/** get the cartesian rank of YLEFT neighbor */
int VCtopology::getYleft_neighbor(){
  return(yleft_neighbor);
}
/** get the cartesian rank of YRIGHT neighbor */
int VCtopology::getYright_neighbor(){
  return(yright_neighbor);
}
/** get the cartesian rank of XLEFT(-) YLEFT(-) SAME Z neighbor  */
int VCtopology::getXleftYleft_neighbor(){
  return(XleftYleft_neighbor);
}
/** cartesian rank of XLEFT(-) YRIGHT(+) SAME Z neighbor */
int VCtopology::getXleftYright_neighbor(){
  return(XleftYright_neighbor);
}
/** cartesian rank of XRIGHT(+) YLEFT(-) SAME Z neighbor */
int VCtopology::getXrightYleft_neighbor(){
  return(XrightYleft_neighbor);
}
/** cartesian rank of XRIGHT(+) YRIGHT(+) SAME Z neighbor */
int VCtopology::getXrightYright_neighbor(){
  return(XrightYright_neighbor);
}

/** get the neighbors in the boundary communicators */
int VCtopology::getLeftNeighbor_COMM_B_LEFT(){
  return(neighborLeft_COMM_B_LEFT);
}
int VCtopology::getRightNeighbor_COMM_B_LEFT(){
  return(neighborRight_COMM_B_LEFT);
}
int VCtopology::getLeftNeighbor_COMM_B_RIGHT(){
  return(neighborLeft_COMM_B_RIGHT);
}
int VCtopology::getRightNeighbor_COMM_B_RIGHT(){
  return(neighborRight_COMM_B_RIGHT);
}
int VCtopology::getLeftNeighbor_COMM_B_BOTTOM(){
  return(neighborLeft_COMM_B_BOTTOM);
}
int VCtopology::getRightNeighbor_COMM_B_BOTTOM(){
  return(neighborRight_COMM_B_BOTTOM);
}
int VCtopology::getLeftNeighbor_COMM_B_TOP(){
  return(neighborLeft_COMM_B_TOP);
}
int VCtopology::getRightNeighbor_COMM_B_TOP(){
  return(neighborRight_COMM_B_TOP);
}

/** if cVERBOSE == true, print to the screen all the comunication */
bool VCtopology::getcVERBOSE(){
   return(cVERBOSE);
}
/** get the coordinates in dir direction of process*/
int VCtopology::getCoordinates(int dir){
   return(coordinates[dir]);
}
/** get Periodicity condition in dir direction */
int VCtopology::getPeriods(int dir){
   return(periods[dir]);
}


// set RefLevelAd                                                                                                                                        
void VCtopology::setRefLevelAdj(int val)
{
  RefLevelAdj= val;
}
// get RefLevelAd                                                                                                                                        
int VCtopology::getRefLevelAdj()
{
  return RefLevelAdj;
}
