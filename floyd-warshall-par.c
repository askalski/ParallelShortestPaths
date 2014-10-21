#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "graph.h"
#include "message_tags.h"
#include "tests.h"
#include "assert.h"

//2do dodać przypadek liczba procesów == 1;

void makeIthRound(int procNum, int thRound, graph ** g, graph ** slots, graph * tempSlot, graph * pivot, int numCells, int cellSize, int numVertices)
{
	int working = 1;
	if(procNum >= numCells) working = 0;

	if(procNum == thRound) //piwotujący
	{
		{
			int i;
			deb("r %d p %d: rozsylam pivot row'a\n", thRound, procNum);

			for(i = 0 ; i < numCells; ++i)
			if(i != procNum)
			{
				//void sendGraph(graph * what, int to, int cellsize);
				sendGraph(g[i], i, cellSize);
			}

			deb("r %d p %d: rozeslalem pivot row'a, zbieram wyniki\n", thRound, procNum);

			for(i = 0 ; i < numCells; ++i)
			if(i != procNum)
			{
				//void broadcastGraph(int who, graph * slot, int cellsize);
				broadcastGraph(i, g[i], cellSize);		
			}
			
			deb("r %d p %d: zebrałem pivot row'a, czekam na pivota+1\n", thRound, procNum);			

			if(thRound < numCells - 1)
				broadcastGraph(thRound+1, pivot, cellSize);
		}
	}
	else // niepiwotujący
	{
		if(working)
		{
			//void receiveGraph(int from, graph * slot, int cellsize);
	
			// pod slots[procNum] zapisujemy swój fragment pivot row'a. tam go aktualizujemy, a potem synchronizujemy.
			receiveGraph(thRound, slots[procNum], cellSize);
		
			deb("r %d p %d: dostaję %d,%d\n", thRound, procNum, slots[procNum]->pos_x, slots[procNum]->pos_y);

			assert(slots[procNum]->pos_y == procNum);
			assert(slots[procNum]->pos_x == thRound);

			//void updatePivotRow(graph * pivot, graph * rowCell, int cellSize);
			//update slot with pivot.

			deb("r %d p %d: licze pivot row'a\n", thRound, procNum);
			updatePivotRow(pivot, slots[procNum], cellSize);
			deb("r %d p %d: policzylem pivot row'a\n", thRound, procNum);
		}

		{
			if(working)
			{
				int i;
				//void broadcastGraph(int who, graph * slot, int cellsize);

				deb("r %d p %d: broadcastuje pivot row'a\n", thRound, procNum);
				for(i = 0 ; i < numCells; ++i)
				if(i != thRound) // dla i = procNum call wysyłania jest taki sam jak call odbierania dla każdego innego.
				{
					broadcastGraph(i, slots[i], cellSize);
				}

				deb("r %d p %d: koniec broadcastu pivot row'a\n", thRound, procNum);
			}
			else
			{
				int i;
				deb("r %d p %d: odbieram niepotrzebne broadcasty\n", thRound, procNum);
				for(i = 0 ; i < numCells; ++i)
				if(i != thRound) // dla i = procNum call wysyłania jest taki sam jak call odbierania dla każdego innego.
				{
					broadcastGraph(i, tempSlot, cellSize);
				}

				deb("r %d p %d: koniec odbioru niepotrzebnych broadcastów\n", thRound, procNum);
			}
		}

		if(working)
		{
				// aktualizuję swój wiersz
			updatePivotCol(pivot, g[thRound], cellSize);
			{
				//void updateGeneral(graph * pivotRowCell, graph * pivotColCell, graph * targetCell, int cellSize)
				int i;
				for(i = 0 ; i < numCells; ++i)
				if(i != thRound)
				{
					updateGeneral(slots[i], g[thRound], g[i], cellSize);
				}			
			}
		}

		if(thRound < numCells-1)
		{
			if(working)
			{
				if(procNum == thRound + 1)
				{
					updatePivot(g[procNum], cellSize);
					broadcastGraph(procNum, g[procNum], cellSize);
				}
				else
				{
					broadcastGraph(thRound+1, pivot, cellSize);
				}
			}
			else
			{
				broadcastGraph(thRound+1, tempSlot, cellSize);
			}
		}

	}
}

int zeroMain(int argc, char * argv[], int procNum, int procSize, int numVertices, int showResults)
{
	graph ** g, **slots, * pivot, * tempSlot;
	double startTime;
	double endTime;

	int numCells;
	int cellSize;
	cellSize = numVertices/procSize + ( numVertices % procSize ? 1 : 0);
	numCells = numVertices/cellSize + ( numVertices % cellSize ? 1 : 0);

	deb("cellSize = %d numCells = %d\n", cellSize, numCells);

	if(numVertices <= 0)
	{
		fprintf(stderr, "Usage: %s [--show-results] <num_vertices>\n", argv[0]);
		MPI_Finalize();
		return 1;
	}

	fprintf(stderr, "Running the Floyd-Warshall algorithm for a graph with %d vertices.\n", numVertices);

	g = (graph**)malloc(sizeof(graph*) * numCells);
	slots = (graph**)malloc(sizeof(graph*) * numCells);

	tempSlot = createGraphSlot(cellSize);
	pivot = createGraphSlot(cellSize);


	if(g == NULL || slots == NULL)
	{
		fprintf(stderr, "Error initializing the graph for the algorithm.\n");
		MPI_Finalize();
		return 2;
	}

	{ // alokacja pamięci pod graf, inicjacja własnego wiersza

		//graph * createGraph(int pos_x, int pos_y, int cellsize, int numVertices);
		//graph * createGraphSlot(int cellsize);
		int i;
		for(i = 0 ; i < numCells; ++i)
			g[i] = createGraphSlot(cellSize);

		for(i = 0 ; i < numCells; ++i)
			slots[i] = createGraphSlot(cellSize);void updatePivotRow(graph * pivot, graph * rowCell, int cellSize);
	}

	//void initGraphParts(graph ** g, int graphSize, int cellSize, int numCells, int actRow);	
	initGraphParts(g, numVertices, cellSize, numCells, 0);

	if(showResults)
	{
		//void printGraphParts(graph ** g, int graphSize, int cellSize, int numCells);
		printGraphParts(g, numVertices, cellSize, numCells);
	}

	{	// wypisywanie debugowe
		int i;
		deb("jestem procesem %d, mam wiersz: \n", procNum);
		for(i = 0 ; i < numCells; ++i)
			deb("[%d, %d] (%d, %d)\n", g[i]->pos_x, g[i]->pos_y, g[i]->size_x, g[i]->size_y);
	}


//	updatePivot(g[0], cellSize);
//	copyGraph(g[0], pivot, cellSize);

	{ // rozsyłam wiersze i ew wypisuje.
		int i, j;
		for(i = 1; i < numCells; ++i)
		{
			initGraphParts(slots, numVertices, cellSize, numCells, i);

			if(showResults)
				printGraphParts(slots, numVertices, cellSize, numCells);

			//void sendGraph(graph * what, int to, int cellsize);
			for(j = 0; j < numCells; ++j)
				sendGraph(slots[j], i, cellSize);

//			sendGraph(pivot, i, cellSize);
		}
	}


	startTime = MPI_Wtime();


	updatePivot(g[0], cellSize);
	copyGraph(g[0], pivot, cellSize);

	{	//rozsyłam zerowego pivota
		int i;
		for(i = 1; i < numCells; ++i)
		{
			sendGraph(pivot, i, cellSize);
		}
	}

	{	//chciałem z tego miejsca pozdrowić mojego kota.
		//void makeIthRound(int procNum, int thRound, graph ** g, graph ** slots, graph * tempSlot, graph * pivot, int numCells, int cellSize, int numVertices)

		int k;
		
		for(k = 0 ; k < numCells; ++k)
		{
			makeIthRound(procNum, k, g, slots, tempSlot, pivot, numCells, cellSize, numVertices);
		}
	}

	deb("p %d skończyłem tańczyć\n", procNum);
  
	endTime = MPI_Wtime();
  
	fprintf(
		stderr,
		"The time required for the Floyd-Warshall algorithm on a %d-node graph with %d process(es): %f.\n",
		numVertices,
		procSize, /* numProcesses */
		endTime - startTime
	);

	deb("zaczynam odbior wynikow\n");
	if(showResults) // agreguje i wywalam wyniki
	{ //void printGraphParts(graph ** g, int graphSize, int cellSize, int numCells);
		int i, j;
		printGraphParts(g, numVertices, cellSize, numCells);
		for(i = 1; i < numCells;++i)
		{
			deb("p %d odbieram wyniki od %d\n", procNum, i);
			for(j = 0 ; j < numCells; ++j)
				receiveGraph(i, slots[j], cellSize);
			deb("p %d odebrałem wyniki od %d\n", procNum, i);
			printGraphParts(slots, numVertices, cellSize, numCells);
			
		}
	}

	{	//zwalnianie pamięci
		int i;

		for(i = 0 ; i < numCells; ++i)
		{
			freeGraph(g[i]);
			freeGraph(slots[i]);
		}
		freeGraph(pivot);
		freeGraph(tempSlot);
	}

	free(g);
	free(slots);

	MPI_Finalize();

	return 0;
}

int nonzeroMain(int argc, char * argv[], int procNum, int procSize, int numVertices, int showResults)
{
	int numCells;
	int cellSize;
	int working;
	graph ** g = NULL, ** slots = NULL, * pivot = NULL, * tempSlot = NULL;
	cellSize = numVertices/procSize + ( numVertices % procSize ? 1 : 0);
	numCells = numVertices/cellSize + ( numVertices % cellSize ? 1 : 0);

	working = 1;
	if(procNum >= numCells) working = 0;

	if(working)
	{
		g = (graph**)malloc(sizeof(graph*) * numCells);
		slots = (graph**)malloc(sizeof(graph*) * numCells);

		pivot = createGraphSlot(cellSize);
	}
	tempSlot = createGraphSlot(cellSize);
	
	if(working)
	{ // alokacja pamięci pod graf

		//graph * createGraph(int pos_x, int pos_y, int cellsize);
		//graph * createGraphSlot(int cellsize);
		int i;
		for(i = 0 ; i < numCells; ++i)
			g[i] = createGraphSlot(cellSize);

		for(i = 0 ; i < numCells; ++i)
			slots[i] = createGraphSlot(cellSize);
	}

	if(working)
	{ // odbieramy swój wiersz.
		int i;
		//void receiveGraph(int from, graph * slot, int cellsize);
		for(i = 0 ; i < numCells; ++i)
		{
			receiveGraph(0, g[i], cellSize);
			//fprintf(stderr, "odbieranie %d porcji\n", i);
		}
		// odbieramy pivota.
		receiveGraph(0, pivot, cellSize);
	}

	if(working)
	{	// wypisywanie debugowe
		int i;
		deb("jestem procesem %d, mam wiersz: \n", procNum);
		for(i = 0 ; i < numCells; ++i)
			deb("[%d, %d] (%d, %d)\n", g[i]->pos_x, g[i]->pos_y, g[i]->size_x, g[i]->size_y);
	}
	else
	{
		deb("jestem procesem %d i nie mam nic do roboty\n", procNum);
	}

	{	//tańczymy
		//void makeIthRound(int procNum, int thRound, graph ** g, graph ** slots, graph * tempSlot, graph * pivot, int numCells, int cellSize, int numVertices)

		int k;
		
		for(k = 0 ; k < numCells; ++k)
		{
			makeIthRound(procNum, k, g, slots, tempSlot, pivot, numCells, cellSize, numVertices);
		}
	}

	deb("p %d skończyłem tańczyć\n", procNum);

	if(working)
	if(showResults) // wysyłam wyniki
	{
		int j;
		for(j = 0 ; j < numCells; ++j)
			sendGraph(g[j], 0, cellSize);
	}

	if(working)
	{	//zwalnianie pamięci
		int i;

		for(i = 0 ; i < numCells; ++i)
		{
			freeGraph(g[i]);
			freeGraph(slots[i]);
		}
		freeGraph(pivot);

		free(g);
		free(slots);
	}
	freeGraph(tempSlot);

	MPI_Finalize();
	return 0;
}

int main(int argc, char * argv[])
{
	int procSize, procNum;
	int numVertices = 0;
	int showResults = 0;

#ifdef USE_RANDOM_GRAPH
#ifdef USE_RANDOM_SEED
	srand(USE_RANDOM_SEED);
#endif
#endif

	MPI_Init(&argc, &argv);

	{
		int i;
		for(i = 1; i < argc; ++i)
		{
			if(strcmp(argv[i], "--show-results") == 0)
			{
				showResults = 1;
			}
			else
			{
				numVertices = atoi(argv[i]);
			}
		}
	}

	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procNum);

  	if(procNum == 0)
		return zeroMain(argc, argv, procNum, procSize, numVertices, showResults);
	else
		return nonzeroMain(argc, argv, procNum, procSize, numVertices, showResults);
}

