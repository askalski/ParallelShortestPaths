#ifndef __graph_h__
#define __graph_h__

#include "stdio.h"

typedef struct 
{
	int size_x, size_y; // opisuje wykorzystany fragment, nie rozmiary tablicy. tablica zawsze jest cellsize x cellsize
	int pos_x, pos_y;
	int * matrix;
} graph;

//graph * createGraph(int pos_x, int pos_y, int cellsize, int numVertices);
graph * createGraphSlot(int cellsize);

void initGraphParts(graph ** g, int graphSize, int cellSize, int numCells, int actRow);
void printGraphParts(graph ** g, int graphSize, int cellSize, int numCells);
void freeGraph(graph * g);

void sendGraph(graph * what, int to, int cellsize);
void receiveGraph(int from, graph * slot, int cellsize);
void broadcastGraph(int who, graph * slot, int cellsize);
size_t getGraphSizeof(int cellsize);

int getCutCellSize(int cellSize, int cellIndex, int numVertices);

void updatePivot(graph * p, int cellSize);
void updatePivotRow(graph * pivot, graph * rowCell, int cellSize);
void updatePivotCol(graph * pivot, graph * colCell, int cellSize);
void updateGeneral(graph * pivotRowCell, graph * pivotColCell, graph * targetCell, int cellSize);
void copyGraph(graph * from, graph * to, int cellSize);

#endif

