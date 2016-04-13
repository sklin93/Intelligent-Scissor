// mex determinePath.cpp to compile
#include "mex.h"
#include <stdlib.h>

//min, max helper function
inline long Min(long a, long b){return a < b ? a : b;}
inline long Max(long a, long b){return a > b ? a : b;}
//helper function to get a 1-D index from 2-D subscripts
inline long getIndex(short xInd, short yInd, short ySize){return long(xInd * ySize + yInd);}

// node structure
struct Node {
    short xPos;
    short yPos;
    long idx;
    float totalCost;
};

// get index of min total cost
long extractMinInd(Node *activeList, long length){
    long    IndOfMin = 0;
    float   curMin   = 10000000000;
    Node  N;
    for (long i = 0; i < length; i++){
        N = *activeList++;
        if (N.totalCost < curMin){
            curMin = N.totalCost;
            IndOfMin = i;   
        }
    }
    return IndOfMin;
}

// Get the Index of the list entry in *activeList idx == indx
long getEqualInd(Node *activeList, long length, long indx){
    Node N;
    for (long i = 0; i < length; i++){
        N = *activeList++;
        if (N.idx == indx) return i;
    }
    return -1;
}

/**********MAIN**********/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    // inputs
    // 1: local costs (a matrix)
    double *localCost = (double*) mxGetData(prhs[0]);
    short xSize = short(mxGetN(prhs[0]));
	short ySize = short(mxGetM(prhs[0]));
    // 2 & 3: xPos & yPos of seed point
    //index in C++ should -1 compared with matlab
    short seedXpos = short(*mxGetPr(prhs[1])) - 1L;//L for not overflow
    short seedYpos = short(*mxGetPr(prhs[2])) - 1L;   
    // 4: radius of pixels from seed point to consider
    double processR;
    if (nrhs < 4) processR = 800; 
    else processR = *mxGetPr(prhs[3]);
    
    const int* imgSize = mxGetDimensions(prhs[0]);
    plhs[0] = mxCreateNumericArray(2, imgSize, mxINT8_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, imgSize, mxINT8_CLASS, mxREAL);
    //pointer to outputs
    char *xPathMap = (char*) mxGetData(plhs[0]);
    char *yPathMap = (char*) mxGetData(plhs[1]);
    
    //variables in use
    long listInd = 0; //for active list length
    short xPosL, xPosU, yPosL, yPosU; //for neighborhood position
    long pixelLeft = Min(long(3.14 * processR * processR + 0.5), long(xSize * ySize));;
    long pixelProcessed = 0;
    float curTotalCost;
    float lengthFactor;
    long Ind;
    long idx;
    
    Node Nq, Nr;
    Node *activeList = (Node*) mxCalloc(10000, sizeof(Node)); 
    char *state = (char*) mxCalloc(long(xSize)*long(ySize), sizeof(char)); //1: processed
    
    // Initialize active list with seed point
    Nq.xPos = seedXpos;
    Nq.yPos = seedYpos;
    Nq.totalCost = 0.0; //seed point cost 0
    Nq.idx = getIndex(seedXpos, seedYpos, ySize);   
    activeList[listInd++] = Nq;

    while (listInd && (pixelProcessed<pixelLeft)){
        pixelProcessed++;
        // find minimal total cost position's index
        Ind = extractMinInd(activeList, listInd);
        Nq = activeList[Ind];
        //remove from active list
        listInd--;
        activeList[Ind] = activeList[listInd];
        //mark state as 'processed'
        state[Nq.idx]  = 1; 
        
        // 3* 3 neighbourhood of Nq
        xPosL = Max(0, Nq.xPos - 1);
        xPosU = Min(xSize - 1, Nq.xPos + 1);
        yPosL = Max(0, Nq.yPos - 1);
        yPosU = Min(ySize - 1, Nq.yPos + 1);
        //loop over neibours
        for (short xPos = xPosL; xPos <= xPosU; xPos++) {
            for (short yPos = yPosL; yPos <= yPosU; yPos++) {
                //this neighbor already processed
                idx = getIndex(xPos, yPos, ySize);
                if (state[idx]) continue;
                // Compute total cost to this neighbour
                if ((abs(xPos - Nq.xPos) + abs(yPos - Nq.yPos)) == 1) 
                    lengthFactor = 1; //horizontal & vertical
                else lengthFactor = 1.41; //diagonal
                curTotalCost = Nq.totalCost + float(localCost[idx])*lengthFactor;         
                Ind = getEqualInd(activeList, listInd, idx);
                if (Ind >= 0){ //in active list
                    Nr = activeList[Ind];
                    // current cost lower
                    if (curTotalCost < Nr.totalCost){
                        Nr.totalCost = curTotalCost;
                        activeList[Ind] = Nr;
                        xPathMap[idx] = char(Nq.xPos - xPos);
                        yPathMap[idx] = char(Nq.yPos - yPos);
                    }
                }else{ //not in active list, place into
                    Nr.xPos = xPos;
                    Nr.yPos = yPos;
                    Nr.idx = idx;
                    Nr.totalCost = curTotalCost;
                    activeList[listInd++] = Nr;
                    xPathMap[idx] = char(Nq.xPos - xPos);
                    yPathMap[idx] = char(Nq.yPos - yPos);
                }
            }
        }
    }
    mxFree(state);
    mxFree(activeList);
}