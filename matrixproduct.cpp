#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <papi.h>
#include <omp.h>

using namespace std;

#define SYSTEMTIME clock_t

// Displays the first 10 elements of the result matrix
void showResultmatrix(double *phc, int n)
{
	int i, j;
	cout << "Result matrix: ";
	for (i = 0; i < 1; i++)
		for (j = 0; j < min(10, n); j++)
			cout << phc[j] << " ";
	cout << endl;
}

// Initializes three matrices...
void initMatrix(double *pha, double *phb, double *phc, int m_ar, int m_br, int m_cr)
{
	int i, j;
	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_ar; j++)
			pha[i * m_ar + j] = (double)1.0;		// ...matrix pha with 1.0

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phb[i * m_br + j] = (double)(i + 1);	// ...matrix phb where phb[i][j] = i + 1

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_br; j++)
			phc[i * m_ar + j] = 0.0;				// ...matrix phc with 0.0
}

// Old call: OnMult(int m_ar, int m_br)
// Added m_cr to represent the number of columns in Matrix A and the number of rows in Matrix B
// Added time to store the execution time of the multiplication
// OnMult -> Algorithm 1
void OnMult(int m_ar, int m_br, int m_cr, double *time) // i -> j -> k
{
	SYSTEMTIME Time1, Time2;

	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	initMatrix(pha, phb, phc, m_ar, m_br, m_cr);

	Time1 = clock();

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_br; j++)
			for (k = 0; k < m_ar; k++)
				phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];

	Time2 = clock();

	*time = (double)(Time2 - Time1) / CLOCKS_PER_SEC;
	cout << "Time: " << *time << " seconds; Size: " << m_ar << ";" << endl;

	showResultmatrix(phc, m_ar * m_ar);

	free(pha);
	free(phb);
	free(phc);
}

// Old call: OnMultLine(int m_ar, int m_br)
// Added m_cr to represent the number of columns in Matrix A and the number of rows in Matrix B
// Added time to store the execution time of the multiplication
// OnMultLine -> Algorithm 2
void OnMultLine(int m_ar, int m_br, int m_cr, double *time) // i -> k -> j
{
    SYSTEMTIME Time1, Time2;

	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	initMatrix(pha, phb, phc, m_ar, m_br, m_cr);

	Time1 = clock();

	for (i = 0; i < m_ar; i++)
		for (k = 0; k < m_ar; k++)
			for (j = 0; j < m_br; j++)
				phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];

	Time2 = clock();

	*time = (double)(Time2 - Time1) / CLOCKS_PER_SEC;
	cout << "Time: " << *time << " seconds; Size: " << m_ar << ";" << endl;

	showResultmatrix(phc, m_ar * m_ar);

	free(pha);
	free(phb);
	free(phc);
    
}

// Old call: OnMultBlock(int m_ar, int m_br, int bkSize)
// Added m_cr to represent the number of columns in Matrix A and the number of rows in Matrix B
// Added time to store the execution time of the multiplication
// OnMultBlock -> Algorithm 3
void OnMultBlock(int m_ar, int m_br, int m_cr, int bkSize, double *time)
{
    SYSTEMTIME Time1, Time2;

	double temp;
	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	initMatrix(pha, phb, phc, m_ar, m_br, m_cr);

	Time1 = clock();

	int numBlocksA = (m_ar + bkSize - 1) / bkSize;
	int numBlocksB = (m_br + bkSize - 1) / bkSize;
	int numBlocksC = (m_cr + bkSize - 1) / bkSize;

	for (int blockA = 0; blockA < numBlocksA; blockA++)
		for (int blockC = 0; blockC < numBlocksC; blockC++)
			for (int blockB = 0; blockB < numBlocksB; blockB++)
			{
				// Verifica que estamos no dentro dos limites da matriz
				int actualBkA = min(bkSize, m_ar - blockA * bkSize);
				int actualBkB = min(bkSize, m_br - blockB * bkSize);
				int actualBkC = min(bkSize, m_cr - blockC * bkSize);

				// Multiplica os blocos
				for (int i = 0; i < actualBkA; i++)
					for (int k = 0; k < actualBkC; k++)
						for (int j = 0; j < actualBkB; j++)

							phc[(i + blockA * bkSize) * m_cr + blockB * bkSize + j] +=
								pha[(blockA * bkSize + i) * m_br + blockC * bkSize + k] *
								phb[(k + blockC * bkSize) * m_cr + blockB * bkSize + j];
			}

	Time2 = clock();

	*time = (double)(Time2 - Time1) / CLOCKS_PER_SEC;
	cout << "Time: " << *time << " seconds; Size: " << m_ar << ";" << endl;

	showResultmatrix(phc, m_ar * m_ar);

	free(pha);
	free(phb);
	free(phc);
}

// OnMultLineParallelV1 -> Algorithm 4
void OnMultLineParallelV1(int m_ar, int m_br, int m_cr, double *time)
{
	SYSTEMTIME Time1, Time2;

	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	initMatrix(pha, phb, phc, m_ar, m_br, m_cr);

	double start_time = omp_get_wtime();

	#pragma omp parallel for private(k, j)
	for (i = 0; i < m_ar; i++)
		for (k = 0; k < m_ar; k++)
			for (j = 0; j < m_br; j++)
				phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];

	double end_time = omp_get_wtime();

	*time = end_time - start_time;
	cout << "Time: " << *time << " seconds; Size: " << m_ar << ";" << endl;

	showResultmatrix(phc, m_ar * m_ar);

	free(pha);
	free(phb);
	free(phc);
}

// OnMultLineParallelV2 -> Algorithm 5
void OnMultLineParallelV2(int m_ar, int m_br, int m_cr, double *time)
{
	SYSTEMTIME Time1, Time2;

	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	initMatrix(pha, phb, phc, m_ar, m_br, m_cr);

	double start_time = omp_get_wtime();

	#pragma omp parallel private(i, k)
	for (i = 0; i < m_ar; i++)
		for (k = 0; k < m_ar; k++)
			#pragma omp for
			for (j = 0; j < m_br; j++)
				phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];

	double end_time = omp_get_wtime();

	*time = end_time - start_time;
	cout << "Time: " << *time << " seconds; Size: " << m_ar << ";" << endl;

	showResultmatrix(phc, m_ar * m_ar);

	free(pha);
	free(phb);
	free(phc);
}

void handle_error (int retval)
{
  printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
  exit(1);
}

void init_papi() {
  int retval = PAPI_library_init(PAPI_VER_CURRENT);
  if (retval != PAPI_VER_CURRENT && retval < 0) {
    printf("PAPI library version mismatch!\n");
    exit(1);
  }
  if (retval < 0) handle_error(retval);

  std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
            << " MINOR: " << PAPI_VERSION_MINOR(retval)
            << " REVISION: " << PAPI_VERSION_REVISION(retval) << "\n";
}


int main (int argc, char *argv[])
{

	if (argc < 6)
	{
		cerr << "Usage: " << argv[0] << " <n> <m> <k> <algorithm> <output_file>\n";
		exit(EXIT_FAILURE);
	}

	int m_ar = atoi(argv[1]);		// n = number of rows of matrix A
	int m_br = atoi(argv[2]); 		// m = number of columns of matrix A (same number of rows of matrix B)
	int m_cr = atoi(argv[3]);		// k = number of columns of matrix B
	int op = atoi(argv[4]);
	string output_file = argv[5];
	int blockSize = (op == 3) ? atoi(argv[6]) : 0;
	double time = 0.0;

	ofstream outfile(output_file, ios::app);
	if (!outfile)
	{
		cerr << "ERROR: cannot open: " << output_file << ".\n";
		exit(EXIT_FAILURE);
	}

	int EventSet = PAPI_NULL;
	long long values[2] = {0, 0};
	int ret;

	ret = PAPI_library_init(PAPI_VER_CURRENT);
	if (ret != PAPI_VER_CURRENT)
		std::cout << "FAIL" << endl;

	ret = PAPI_create_eventset(&EventSet);
	if (ret != PAPI_OK)
		cout << "ERROR: create eventset" << endl;

	ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
	if (ret != PAPI_OK)
		cout << "ERROR: PAPI_L1_DCM" << endl;

	ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
	if (ret != PAPI_OK)
		cout << "ERROR: PAPI_L2_DCM" << endl;

	ret = PAPI_start(EventSet);
	if (ret != PAPI_OK)
		cout << "ERROR: Start PAPI" << endl;

	switch (op)
	{
	case 1:
		OnMult(m_ar, m_br, m_cr, &time);
		break;
	case 2:
		OnMultLine(m_ar, m_br, m_cr, &time);
		break;
	case 3:
		OnMultBlock(m_ar, m_br, m_cr, blockSize, &time);
		break;
	case 4:
		OnMultLineParallelV1(m_ar, m_br, m_cr, &time);
		break;
	case 5:
		OnMultLineParallelV2(m_ar, m_br, m_cr, &time);
		break;
	default:
		cerr << "Usage: " << argv[0] << " <n> <m> <k> <algorithm> <output_file>\n";
		exit(EXIT_FAILURE);
	}

	ret = PAPI_stop(EventSet, values);
	if (ret != PAPI_OK)
		cout << "ERROR: Stop PAPI" << endl;

	// Store the results
	outfile << op << ","
			<< m_ar << ","
			<< blockSize << ","
			<< time << ","
			<< values[0] << ","
			<< values[1] << endl;

	// TODO: I removed the PAPI cleanup

}