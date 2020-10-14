#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <windows.h>
#include <chrono>
#include <omp.h>
#include <thread>
#include <sysinfoapi.h>
#include <intrin.h>
#include <iostream>
#include <fstream>
#include <ctime>    


//namespace improved {
int number_bacteria;
char** bacteria_name;
double** comparisons; // results of CompareBacteria
int NUM_THREADS; // num threads for OMP to use
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

void Init()
{
	M2 = 1;
	for (int i = 0; i < LEN - 2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
		M2 *= AA_NUMBER;
	M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
	M = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
}

class Bacteria
{
private:
	long* vector;
	long* second;
	long one_l[AA_NUMBER];
	long indexs;
	long total;
	long total_l;
	long complement;

	void InitVectors()
	{
		vector = new long[M];
		second = new long[M1];
		memset(vector, 0, M * sizeof(long));
		memset(second, 0, M1 * sizeof(long));
		memset(one_l, 0, AA_NUMBER * sizeof(long));
		total = 0;
		total_l = 0;
		complement = 0;
	}

	void init_buffer(char* buffer)
	{
		complement++;
		indexs = 0;
		for (int i = 0; i < LEN - 1; i++)
		{
			short enc = encode(buffer[i]);
			one_l[enc]++;
			total_l++;
			indexs = indexs * AA_NUMBER + enc;
		}
		second[indexs]++;
	}

	void cont_buffer(char ch)
	{
		short enc = encode(ch);
		one_l[enc]++;
		total_l++;
		long index = indexs * AA_NUMBER + enc;
		vector[index]++;
		total++;
		indexs = (indexs % M2) * AA_NUMBER + enc;
		second[indexs]++;
	}

public:
	long count;
	double* tv;
	long* ti;

	Bacteria(char* filename)
	{
		FILE* bacteria_file;
		errno_t OK = fopen_s(&bacteria_file, filename, "r");

		if (OK != 0)
		{
			fprintf(stderr, "Error: failed to open file %s\n", filename);
			exit(1);
		}

		InitVectors();

		char ch;
		while ((ch = fgetc(bacteria_file)) != EOF)
		{
			if (ch == '>')
			{
				while (fgetc(bacteria_file) != '\n'); // skip rest of line

				char buffer[LEN - 1];
				fread(buffer, sizeof(char), LEN - 1, bacteria_file);
				init_buffer(buffer);
			}
			else if (ch != '\n')
				cont_buffer(ch);
		}

		long total_plus_complement = total + complement;
		double total_div_2 = total * 0.5;
		int i_mod_aa_number = 0;
		int i_div_aa_number = 0;
		long i_mod_M1 = 0;
		long i_div_M1 = 0;

		double one_l_div_total[AA_NUMBER];
		for (int i = 0; i < AA_NUMBER; i++)
			one_l_div_total[i] = (double)one_l[i] / total_l;

		double* second_div_total = new double[M1];
		for (int i = 0; i < M1; i++)
			second_div_total[i] = (double)second[i] / total_plus_complement;

		count = 0;
		double* t = new double[M];

		for (long i = 0; i < M; i++)
		{
			double p1 = second_div_total[i_div_aa_number];
			double p2 = one_l_div_total[i_mod_aa_number];
			double p3 = second_div_total[i_mod_M1];
			double p4 = one_l_div_total[i_div_M1];
			double stochastic = (p1 * p2 + p3 * p4) * total_div_2;

			if (i_mod_aa_number == AA_NUMBER - 1)
			{
				i_mod_aa_number = 0;
				i_div_aa_number++;
			}
			else
				i_mod_aa_number++;

			if (i_mod_M1 == M1 - 1)
			{
				i_mod_M1 = 0;
				i_div_M1++;
			}
			else
				i_mod_M1++;

			if (stochastic > EPSILON)
			{
				t[i] = (vector[i] - stochastic) / stochastic;
				count++;
			}
			else
				t[i] = 0;
		}

		delete second_div_total;
		delete vector;
		delete second;

		tv = new double[count];
		ti = new long[count];

		int pos = 0;
		for (long i = 0; i < M; i++)
		{
			if (t[i] != 0)
			{
				tv[pos] = t[i];
				ti[pos] = i;
				pos++;
			}
		}
		delete t;

		fclose(bacteria_file);
	}
};

void ReadInputFile(const char* input_name)
{
	FILE* input_file;
	errno_t OK = fopen_s(&input_file, input_name, "r");

	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
		exit(1);
	}

	fscanf_s(input_file, "%d", &number_bacteria);
	bacteria_name = new char* [number_bacteria];

	for (long i = 0;i < number_bacteria;i++)
	{
		char name[10];
		fscanf_s(input_file, "%s", name, 10);
		bacteria_name[i] = new char[20];
		sprintf_s(bacteria_name[i], 20, "data/%s.faa", name);
	}
	fclose(input_file);
}

double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
	double correlation = 0;
	double vector_len1 = 0;
	double vector_len2 = 0;
	long p1 = 0;
	long p2 = 0;

	while (p1 < b1->count && p2 < b2->count)
	{
		long n1, n2;
		n1 = b1->ti[p1];
		n2 = b2->ti[p2];

		if (n1 < n2)
		{
			double t1;
			t1 = b1->tv[p1];
			vector_len1 += (t1 * t1);
			p1++;
		}
		else if (n2 < n1)
		{
			double t2;
			t2 = b2->tv[p2];
			p2++;
			vector_len2 += (t2 * t2);
		}
		else
		{
			double t1, t2;
			t1 = b1->tv[p1++];
			t2 = b2->tv[p2++];
			vector_len1 += (t1 * t1);
			vector_len2 += (t2 * t2);
			correlation += t1 * t2;
		}
	}

	while (p1 < b1->count)
	{
		long n1 = b1->ti[p1];
		double t1 = b1->tv[p1++];
		vector_len1 += (t1 * t1);
	}
	while (p2 < b2->count)
	{
		long n2 = b2->ti[p2];
		double t2 = b2->tv[p2++];
		vector_len2 += (t2 * t2);
	}

	return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

struct Comparison {
	int bacteria1;
	int bacteria2;
	double correlation;
};

void CompareAllBacteria()
{
	Bacteria** b = new Bacteria * [number_bacteria];

	// Init comparisons[][] array
	comparisons = new double* [number_bacteria];
	for (int i = 0; i < number_bacteria; ++i)
		comparisons[i] = new double[number_bacteria];

	int progress = 1;
#pragma omp parallel for schedule(dynamic,1) num_threads(NUM_THREADS)
	for (int i = 0; i < number_bacteria; i++)
	{
		b[i] = new Bacteria(bacteria_name[i]);
		
#pragma omp critical
		printf("loaded %d of %d on thread %d\r", progress, number_bacteria, omp_get_thread_num());
#pragma omp atomic
		progress++;
	
	}

	printf("\nRunning comparsions...\n");

	progress = 0;
#pragma omp parallel for schedule(dynamic,1) num_threads(NUM_THREADS)
	for (int i = 0; i < number_bacteria - 1; i++) {
		for (int j = i + 1; j < number_bacteria; j++)
		{
			comparisons[i][j] = CompareBacteria(b[i], b[j]);
		}
#pragma omp critical
		printf("computed %d of %d\r", progress, number_bacteria);
#pragma omp atomic
		progress++;
	}


	// print output in order
	for (int i = 0; i < number_bacteria - 1; i++)
		for (int j = i + 1; j < number_bacteria; j++)
		{
			printf("%2d %2d -> %.20lf\n", i, j, comparisons[i][j]);
		}
}

/**
* Get number of threads to use from argv and set to global var consumed by omp pragma
  If no arg is defined then will try and use max cores for machine. 
  If for some reason the program fails to get the machine core count, default to 2.
*/
void setNumThreadsFromArgs(int argc, char** argv) {
	if (argc == 1) {
		const auto processor_count = std::thread::hardware_concurrency();
		if (processor_count == 0) {
			NUM_THREADS = 2;
			printf("Using MIN (%d) threads\n", NUM_THREADS);
		}
		NUM_THREADS = processor_count;
		printf("Using MAX (%d) threads\n", NUM_THREADS);
	}
	else {
		if (strcmp(argv[1], "--num-threads") == 0) {
			int first_arg = (int)strtol(argv[2], NULL, 10);
			NUM_THREADS = first_arg;
			printf("Using USER SPEC (%d) threads\n", NUM_THREADS);
		}
		else {
			printf("Usage: CvTree.exe --num-threads [number]");
			exit(0);
		}
		
	}
}

/**
* Write the results of CompareBacteria to a file using provided filename
*/
void writeOutput(std::string filename) {
	printf("Writing output to file...\n");

#pragma warning(suppress : 4996) // ignore compiler warning C4996 'may be unsafe'
	FILE* pFile = fopen(filename.c_str(), "w");

	for (int i = 0; i < number_bacteria - 1; i++) {
		for (int j = i + 1; j < number_bacteria; j++) {
			fprintf(pFile, "%2d %2d -> %.20lf\n", i, j, comparisons[i][j]);
		}
	}
	fclose(pFile);

	printf("Wrote output to %s\n", filename.c_str());
}

// Writes num cores and execution time to a file using date and time as filename
void writeTestResult(std::chrono::duration<double> *time) {
	int CPUInfo[4] = { -1 };
	unsigned   nExIds, i = 0;
	char CPUBrandString[0x40];
	// Get the information associated with each extended ID.
	__cpuid(CPUInfo, 0x80000000);
	nExIds = CPUInfo[0];
	for (i = 0x80000000; i <= nExIds; ++i)
	{
		__cpuid(CPUInfo, i);
		// Interpret CPU brand string
		if (i == 0x80000002)
			memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
		else if (i == 0x80000003)
			memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
		else if (i == 0x80000004)
			memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
	}
	//string includes manufacturer, model and clockspeed
	SYSTEM_INFO sysInfo;
	GetSystemInfo(&sysInfo);

	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);

	printf("\nWriting results file...\n");

	auto time_now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

	time_t t = std::time(0);   // get time now
#pragma warning(suppress : 4996) // ignore compiler warning C4996 'may be unsafe'
	struct tm* now = localtime(&t);

	char buffer[80];
	strftime(buffer, 80, "%Y-%m-%d-%H-%M-%S.txt", now);
	std::string filename = buffer;

	std::ofstream file;

	file.open(filename);
	file << "threads used" << NUM_THREADS << std::endl;
	file << "time elapsed " << time->count() << " seconds." << std::endl;	
	file << "CPU Type: " << CPUBrandString << std::endl;
	file << "Number of Cores: " << sysInfo.dwNumberOfProcessors << std::endl;
	file << "Total System Memory: " << (statex.ullTotalPhys / 1024) / 1024 << "MB" << std::endl;


	file.close();

	printf("Wrote output to %s\n", filename.c_str());
}

int main(int argc, char* argv[])
{
	using namespace std::chrono;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
	setNumThreadsFromArgs(argc, argv);
	Init();
	ReadInputFile("list.txt");
	CompareAllBacteria();

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "time elapsed " << time_span.count() << " seconds.";

	//writeOutput("outPutParallel.txt");
	writeTestResult(&time_span);
	return 0;
}
//}