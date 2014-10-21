#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include "assert.h"
#include "string.h"
#include "mpi.h"
#include "message_tags.h"
#include "tests.h"

// i - wiersz, j - kolumna

inline size_t getGraphSizeof(int cellsize)
{
	return sizeof(graph) + sizeof(int) * cellsize * cellsize;
}

graph * createGraphSlot(int cellsize)
{
	graph * g = (graph *)malloc(getGraphSizeof(cellsize));
	if(g == NULL)
		return NULL;

	g->matrix = (int*)((void*)g + sizeof(graph));

	return g;
}

void freeGraph(graph * g)
{
	free(g);
}


void initGraphParts(graph ** g, int graphSize, int cellSize, int numCells, int actRow)
{
	int i,j;
	int * gr;

	for(i = 0 ; i < numCells; ++i)
	{
		g[i]->pos_x = actRow;
		g[i]->pos_y = i;

		//int getCutCellSize(int cellSize, int cellIndex, int numVertices)
		g[i]->size_x = getCutCellSize(cellSize, actRow, graphSize);
		g[i]->size_y = getCutCellSize(cellSize, i, graphSize);
		
		//fprintf(stderr, "komórka [%d,%d] rozmiarów [%d, %d]\n", g[j]->pos_x, g[j]->pos_y, g[j]->size_x, g[j]->size_y);
	}

	for(i = actRow*cellSize ; i < (actRow+1)*cellSize; ++i)
	{
		for(j = 0; j < graphSize; ++j)
		{
			gr = g[j/cellSize]->matrix;
			gr[(i % cellSize)*cellSize + (j%cellSize)] = i == j ? 0 :
#ifndef USE_RANDOM_GRAPH
				(i - j == 1 || j - i == 1 ? 1 : n + 5);
#else
				(rand() & 8191) + 1;
#endif
		}
	}
}

void printGraphParts(graph ** g, int graphSize, int cellSize, int numCells)
{
	int i,j;
	int * gr;
	for(i = 0 ; i < g[0]->size_x; ++i)
	{
		gr = g[0]->matrix;
		printf("%4d", gr[(i % cellSize)*cellSize + 0]);
		for(j = 1; j < graphSize; ++j)
		{
			gr = g[j/cellSize]->matrix;
			printf(" %4d", gr[(i % cellSize)*cellSize + (j%cellSize)]);
		}
		printf("\n");
	}
}

void sendGraph(graph * what, int to, int cellsize)
{
	MPI_Send(
		what,
		getGraphSizeof(cellsize),
		MPI_BYTE,
		to,
		FWP_MATRIX_TAG,
		MPI_COMM_WORLD
	);
}

void receiveGraph(int from, graph * slot, int cellsize)
{
	MPI_Status status;

	MPI_Recv(
		(void*)slot,
		getGraphSizeof(cellsize),
		MPI_BYTE,
		from,
		FWP_MATRIX_TAG,
		MPI_COMM_WORLD,
		&status);

	slot->matrix = (int*)((void*)slot + sizeof(graph));
}

void broadcastGraph(int who, graph * slot, int cellsize)
{
	MPI_Bcast(
		(void*)slot,
		getGraphSizeof(cellsize),
		MPI_BYTE,
		who,
		MPI_COMM_WORLD);	

	slot->matrix = (int*)((void*)slot + sizeof(graph));
}

int getCutCellSize(int cellSize, int cellIndex, int numVertices)
{
	if((cellIndex+1)*cellSize > numVertices) return numVertices % cellSize; // double checked!
	else return cellSize;
}

void updatePivot(graph * p, int cellSize)
{
	int k, i, j;
	int * m = p->matrix;
	int n = p->size_x;
	assert(p->size_x == p->size_y); //pivot!

	deb("pivot-aktualizuje [%d, %d] (%d, %d)\n", p->pos_x, p->pos_y, p->size_x, p->size_y);

	for (k = 0; k < n; ++k)
	{
		for (i = 0; i < n; ++i)
		{
			for (j = 0; j < n; ++j)
			{
				int pathSum = m[i*cellSize + k] + m[k*cellSize + j];
				if (m[i*cellSize + j] > pathSum)
					m[i*cellSize + j] = pathSum;
			}
		}
	}
}

void updatePivotRow(graph * pivot, graph * rowCell, int cellSize)
{
	int k, i, j;

	assert(pivot->pos_x == rowCell->pos_x);
	assert(rowCell->size_x == pivot->size_x);
	assert(pivot->size_x == pivot->size_y);

	deb("pivot-row-aktualizuje [%d, %d] (%d, %d) pivotem [%d, %d] (%d, %d)\n",
			rowCell->pos_x, rowCell->pos_y, rowCell->size_x, rowCell->size_y,
			pivot->pos_x, pivot->pos_y, pivot->size_x, pivot->size_y);

	for(k = 0; k < pivot->size_x; ++k) //pivot sx = sy
	{
		for(i = 0; i < rowCell->size_x; ++i)
		{
			for(j = 0 ; j < rowCell->size_y; ++j)
			{
				int pathSum = pivot->matrix[i*cellSize + k] + rowCell->matrix[k*cellSize + j];
				if(rowCell->matrix[i*cellSize + j] > pathSum)
				   rowCell->matrix[i*cellSize + j] = pathSum;
			}
		}
	}
}

void updatePivotCol(graph * pivot, graph * colCell, int cellSize)
{
	int k, i, j;

	assert(pivot->pos_y == colCell->pos_y);
	assert(colCell->size_y == pivot->size_y);
	assert(pivot->size_x == pivot->size_y);

	deb("pivot-col-aktualizuje [%d, %d] (%d, %d) pivotem [%d, %d] (%d, %d)\n",
			colCell->pos_x, colCell->pos_y, colCell->size_x, colCell->size_y,
			pivot->pos_x, pivot->pos_y, pivot->size_x, pivot->size_y);

	for(k = 0; k < pivot->size_x; ++k)
	{
		for(i = 0; i < colCell->size_x; ++i)
		{
			for(j = 0 ; j < colCell->size_y; ++j)
			{
				int pathSum = colCell->matrix[i*cellSize + k] + pivot->matrix[k*cellSize + j];
				if(colCell->matrix[i*cellSize + j] > pathSum)
				   colCell->matrix[i*cellSize + j] = pathSum;
			}
		}
	}
}


void updateGeneral(graph * pivotRowCell, graph * pivotColCell, graph * targetCell, int cellSize)
{
	int k, i, j, limk;

	deb("aktualizuje targetCell [%d, %d] rozmiarów (%d, %d) za pomocą pRC [%d, %d] (%d, %d) i pCC [%d, %d] (%d, %d)\n",
		targetCell->pos_x, targetCell->pos_y, targetCell->size_x, targetCell->size_y,
		pivotRowCell->pos_x, pivotRowCell->pos_y, pivotRowCell->size_x, pivotRowCell->size_y,
		pivotColCell->pos_x, pivotColCell->pos_y, pivotColCell->size_x, pivotColCell->size_y);


	assert(pivotRowCell->size_y == targetCell->size_y);
	assert(pivotColCell->size_x == targetCell->size_x);

	//gdzieś tu czai się błąd.
	
	//limk = pivotRowCell->size_x;

	limk = pivotColCell->size_y;
	if(limk >= pivotRowCell->size_x) limk = pivotRowCell->size_x;

	for(k = 0; k < limk; ++k)
	{
		for(i = 0; i < pivotColCell->size_x; ++i)
		{
			for(j = 0 ; j < pivotRowCell->size_y; ++j)
			{
				//postaw warunek na zakres tablicy!!!

				//if(k >= pivotColCell->size_y) continue;
				//if(k >= pivotRowCell->size_x) continue;

				int pathSum = pivotColCell->matrix[i*cellSize + k] + pivotRowCell->matrix[k*cellSize + j];
				if(targetCell->matrix[i*cellSize + j] > pathSum)
				   targetCell->matrix[i*cellSize + j] = pathSum;
			}
		}
	}
}

void copyGraph(graph * from, graph * to, int cellSize)
{
	to->pos_x = from->pos_x;
	to->pos_y = from->pos_y;
	to->size_x = from->size_x;
	to->size_y = from->size_y;

	memcpy((void*)to->matrix, (void*)from->matrix, sizeof(int) * cellSize * cellSize);
}


