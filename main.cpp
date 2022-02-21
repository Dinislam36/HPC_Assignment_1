#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <CL/cl.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <cpuid.h>
#include <iostream>
#include <unistd.h>
#include <sysinfoapi.h>
#include <windows.h>

#define MAX_VAL 2;

// Random
double next()
{
    return ((double)rand() / (double)RAND_MAX);
}

// Task 3
struct regMatrix{
    int N;
    double** matrix;
};

void AllocateRegMatrix(int N, regMatrix &mtx){
    mtx.N = N;
    mtx.matrix = (double**)malloc(N * sizeof(double*));
    for(int i = 0; i < N; i++){
        mtx.matrix[i] = (double*)malloc(N * sizeof(double));
    }
}

void FreeRegMatrix(regMatrix &mtx){
    for(int i = 0; i < mtx.N; i++){
        std::free(mtx.matrix[i]);
    }
    std::free(mtx.matrix);
}

// Task 4
struct crsMatrix{
    // Matrix size (N x N)
    int N;
    // Number of non-zero elements
    int NZ;
    // Array of values (size NZ)
    double* Value;
    // Array of column numbers (size NZ)
    int* Col;
    // Array of row indexes (size N + 1)
    int* row_index;
};

void AllocateCRSMatrix(int N, int NZ, crsMatrix &mtx){
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.Value = (double*)malloc(sizeof(double) * NZ);
    mtx.Col = (int*)malloc(sizeof(int) * NZ);
    mtx.row_index = (int*)malloc(sizeof(int) * (N+1));
}

void FreeCRSMatrix(crsMatrix &mtx){
    std::free(mtx.Value);
    std::free(mtx.Col);
    std::free(mtx.row_index);
}

// Task 5
void AllocateVector(int N, double **vec){
    *vec = (double*)malloc(sizeof(double) * N);
}

void FreeVector(double **vec){
    std::free(*vec);
}

// Task 6
void GenerateRegularCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
    int i, j, k, f, tmp, notNull, c;

    srand(seed);

    notNull = cntInRow * N;
    AllocateCRSMatrix(N, notNull, mtx);

    for(i = 0; i < N; i++)
    {
        for(j = 0; j < cntInRow; j++)
        {
            do
            {
                mtx.Col[i * cntInRow + j] = rand() % N;
                f = 0;
                for (k = 0; k < j; k++)
                    if (mtx.Col[i * cntInRow + j] == mtx.Col[i * cntInRow + k])
                        f = 1;
            } while (f == 1);
        }
        for (j = 0; j < cntInRow - 1; j++)
            for (k = 0; k < cntInRow - 1; k++)
                if (mtx.Col[i * cntInRow + k] > mtx.Col[i * cntInRow + k + 1])
                {
                    tmp = mtx.Col[i * cntInRow + k];
                    mtx.Col[i * cntInRow + k] = mtx.Col[i * cntInRow + k + 1];
                    mtx.Col[i * cntInRow + k + 1] = tmp;
                }
    }

    for (i = 0; i < cntInRow * N; i++)
        mtx.Value[i] = next() * MAX_VAL;

    c = 0;
    for (i = 0; i <= N; i++)
    {
        mtx.row_index[i] = c;
        c += cntInRow;
    }
}

void GenerateSpecialCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
    srand(seed);
    double end = pow((double)cntInRow, 1.0 / 3.0);
    double step = end / N;

    std::vector<int> *columns = new std::vector<int>[N];
    int NZ = 0;

    for (int i = 0; i < N; i++)
    {
        int rowNZ = int(pow((double(i + 1) * step), 3) + 1);
        NZ += rowNZ;
        int num1 = (rowNZ - 1) / 2;
        int num2 = rowNZ - 1 - num1;

        if (rowNZ != 0)
        {
            if (i < num1)
            {
                num2 += num1 - i;
                num1 = i;
                for(int j = 0; j < i; j++)
                    columns[i].push_back(j);
                columns[i].push_back(i);
                for(int j = 0; j < num2; j++)
                    columns[i].push_back(i + 1 + j);
            }
            else
            {
                if (N - i - 1 < num2)
                {
                    num1 += num2 - (N - 1 - i);
                    num2 = N - i - 1;
                }
                for (int j = 0; j < num1; j++)
                    columns[i].push_back(i - num1 + j);
                columns[i].push_back(i);
                for (int j = 0; j < num2; j++)
                    columns[i].push_back(i + j + 1);
            }
        }
    }

    AllocateCRSMatrix(N, NZ, mtx);

    int count = 0;
    int sum = 0;
    for (int i = 0; i < N; i++)
    {
        mtx.row_index[i] = sum;
        sum += columns[i].size();
        for (unsigned int j = 0; j < columns[i].size(); j++)
        {
            mtx.Col[count] = columns[i][j];
            mtx.Value[count] = next();
            count++;
        }
    }
    mtx.row_index[N] = sum;

    delete [] columns;
}

void GenerateRegularMatrix(int seed, int N, regMatrix& mtx){
    int i, j;

    srand(seed);

    AllocateRegMatrix(N, mtx);

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            mtx.matrix[i][j] = next() * MAX_VAL;
        }
    }
}

// Task 7
void CRSToRegular(crsMatrix crs, regMatrix& reg){
    AllocateRegMatrix(crs.N, reg);
    for(int i = 0; i < crs.N; i++){
        for(int j = 0; j < crs.N; j++){
            reg.matrix[i][j] = 0;
        }
    }

    int row_index[crs.N+1];
    for(int i = 0; i < crs.N; i++){
        row_index[i] = crs.row_index[i];
    }
    row_index[crs.N] = crs.NZ;
    for(int i = 0; i < crs.N; i++){
        for(int j = crs.row_index[i]; j < crs.row_index[i+1]; j++){
            reg.matrix[i][crs.Col[j]] = crs.Value[j];
        }
    }


}

// Task 8
void SeqRegMult(regMatrix A, double *x, double *b, double &time){
    double val;
    int i, j;
    double start = omp_get_wtime();

    for(i = 0; i < A.N; i++){
        val = 0;
        for(j = 0; j < A.N; j++){
            val = val + A.matrix[i][j] * x[j];
        }
        b[i] = val;
    }
    time = omp_get_wtime() - start;
}

// Task 9
void SeqCRSMult(crsMatrix A, double *x, double *b, double &time){
    double sum;

    for(int i = 0; i < A.N; i++){
        sum = 0;
        for(int j = A.row_index[i]; j < A.row_index[i + 1]; j++){
            sum = sum + A.Value[j] * x[A.Col[j]];
        }
        b[i] = sum;
    }
}

// Task 10
int CompareVectors(double* vec1, double* vec2, int n, double &diff)
{
    diff = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (diff < fabs(vec1[i] - vec2[i]))
        {
            diff = fabs(vec1[i] - vec2[i]);
        }
    }
    return 0;
}

int CompareRegMatrix(regMatrix A, regMatrix B, double &diff){
    if(A.N != B.N)
        return 1;
    diff = 0.0;
    int n = A.N;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(diff < fabs(A.matrix[i][j] - B.matrix[i][j])){
                diff = fabs(A.matrix[i][j] - B.matrix[i][j]);
            }
        }
    }
    return 0;
}

int CompareCRSMatrix(crsMatrix A, crsMatrix B, double &diff){
    if(A.N != B.N or A.NZ != B.NZ)
        return 1;
    diff = 0.0;
    int nz = A.NZ;
    for(int i = 0; i < nz; i++){
        if(diff < fabs(A.Value[i] - B.Value[i])){
            diff = fabs(A.Value[i] - B.Value[i]);
        }
    }
    return 0;
}

// Task 12
void OmpRegMult(regMatrix A, double *x, double *b, int numThreads, double &time){
    int i, j;
    double val;

    double start = omp_get_wtime();

    #pragma omp parallel shared(A, x, b) private(val, j) num_threads(numThreads)
    {
    #pragma omp for
        for(i = 0; i < A.N; i++){
            val = 0;
            for(j = 0; j < A.N; j++){
                val += A.matrix[i][j] * x[j];
            }
            b[i] = val;
        }
    }
    time = omp_get_wtime() - start;
}

// Task 14
void OmpCRSMult(crsMatrix A, double *x, double *b, int numThreads, double &time){
    double sum;
    int i, j;

    double start = omp_get_wtime();
    #pragma omp parallel shared(A, x, b) private(sum, j) num_threads(numThreads)
    {
    #pragma omp for
        for (i = 0; i < A.N; i++) {
            sum = 0;
            for (j = A.row_index[i]; j < A.row_index[i + 1]; j++) {
                sum = sum + A.Value[j] * x[A.Col[j]];
            }
            b[i] = sum;
        }
    }
    time = omp_get_wtime() - start;
}

// Task 16
int SeqRegMatMult(regMatrix A, regMatrix B, regMatrix &C, double &time){
    if(A.N != B.N)
        return 1;

    int N = A.N;

    AllocateRegMatrix(N, C);
    double start = omp_get_wtime();
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            C.matrix[i][j] = 0;
            for(int k = 0; k < N; k++){
                C.matrix[i][j] += A.matrix[i][k] * B.matrix[k][j];
            }
        }
    }
    time = omp_get_wtime() - start;
    return 0;
}

double Transpose2(crsMatrix imtx, crsMatrix &omtx)
{
    int i, j;

    double start = omp_get_wtime();

    AllocateCRSMatrix(imtx.N, imtx.NZ, omtx);

    memset(omtx.row_index, 0, (imtx.N + 1) * sizeof(int));
    for (i = 0; i < imtx.NZ; i++)
        omtx.row_index[imtx.Col[i] + 1]++;

    int S = 0;
    for (i = 1; i <= imtx.N; i++)
    {
        int tmp = omtx.row_index[i];
        omtx.row_index[i] = S;
        S = S + tmp;
    }

    for (i = 0; i < imtx.N; i++)
    {
        int j1 = imtx.row_index[i];
        int j2 = imtx.row_index[i+1];
        int Col = i; // Столбец в AT - строка в А
        for (j = j1; j < j2; j++)
        {
            double V = imtx.Value[j];  // Значение
            int RIndex = imtx.Col[j];  // Строка в AT
            int IIndex = omtx.row_index[RIndex + 1];
            omtx.Value[IIndex] = V;
            omtx.Col  [IIndex] = Col;
            omtx.row_index[RIndex + 1]++;
        }
    }


    return omp_get_wtime() - start;
}

// Task 17
int SeqCRSMatMult(crsMatrix A, crsMatrix B_1, crsMatrix &C, double &time){
    crsMatrix B;
    Transpose2(B_1, B);
    if (A.N != B.N) return 1;
    int N = A.N;
    std::vector<int> columns;
    std::vector<double> values;
    std::vector<int> row_index;
    double start = omp_get_wtime();
    int rowNZ; row_index.push_back(0);
    for (int i = 0; i < N; i++) {
        rowNZ = 0;
        for (int j = 0; j < N; j++) {
            double sum = 0;
// dot product of the rows А и BT
            for (int k = A.row_index[i]; k < A.row_index[i + 1]; k++)
                for (int l = B.row_index[j]; l < B.row_index[j + 1]; l++)
                    if (A.Col[k] == B.Col[l])
                    {
                        sum += A.Value[k] * B.Value[l];
                        break;
                    }
            if (fabs(sum) > 0.00005) {
                columns.push_back(j); values.push_back(sum); rowNZ++;
            }
        }
        row_index.push_back(rowNZ + row_index[i]);
    }
    AllocateCRSMatrix(N, columns.size(), C);
    for (int j = 0; j < columns.size(); j++) {
        C.Col[j] = columns[j]; C.Value[j] = values[j];
    }
    for(int i = 0; i <= N; i++)
        C.row_index[i] = row_index[i];
    clock_t finish = clock();
    time = omp_get_wtime() - start;
    return 0;

}

// Task 19
int OmpRegMatMult(regMatrix A, regMatrix B, regMatrix &C, int numThreads, double &time){
    if(A.N != B.N)
        return 1;

    int N = A.N;

    AllocateRegMatrix(N, C);
    double start = omp_get_wtime();
    int i, j, k;

    #pragma omp parallel shared(A, B, C) private(j, k) num_threads(numThreads)
    {
    #pragma omp for
        for(i = 0; i < N; i++) {
            for(j = 0; j < N; j++) {
                C.matrix[i][j] = 0;
                for(k = 0; k < N; k++) {
                    C.matrix[i][j] += A.matrix[i][k] * B.matrix[k][j];
                }
            }
        }
    }
    time = omp_get_wtime() - start;
    return 0;
}

// Task 21
int OmpCRSMatMult(crsMatrix A, crsMatrix B_1, crsMatrix &C, int numThreads, double &time){
    crsMatrix B;
    Transpose2(B_1, B);
    if (A.N != B.N) return 1;
    int N = A.N;
    std::vector<int> columns;
    std::vector<double> values;
    std::vector<int> row_index;
    double start = omp_get_wtime();
    int rowNZ; row_index.push_back(0);

    int i, j, k, l;
    double sum;

        for (i = 0; i < N; i++) {
            rowNZ = 0;
            #pragma omp parallel shared(A, B, rowNZ, columns, values) private(sum, k, l) num_threads(numThreads)
            {
                #pragma omp for
                for (j = 0; j < N; j++) {
                    sum = 0;
                    // dot product of the rows А и BT
                    for (k = A.row_index[i]; k < A.row_index[i + 1]; k++)
                        for (l = B.row_index[j]; l < B.row_index[j + 1]; l++)
                            if (A.Col[k] == B.Col[l]) {
                                sum += A.Value[k] * B.Value[l];
                                break;
                            }
                    #pragma omp critical(setSum)
                    {
                        if (fabs(sum) > 0.00005) {
                            columns.push_back(j);
                            values.push_back(sum);
                            rowNZ++;
                        }
                    }
                }
            }
            row_index.push_back(rowNZ + row_index[i]);
    }

    AllocateCRSMatrix(N, columns.size(), C);

    #pragma omp parallel shared(C, columns, values) num_threads(numThreads)
    {
        #pragma omp for
        for (j = 0; j < columns.size(); j++) {
            C.Col[j] = columns[j];
            C.Value[j] = values[j];
        }
    }

    #pragma omp parallel shared(C, row_index) num_threads(numThreads)
    {
        #pragma omp for
        for (i = 0; i <= N; i++)
            C.row_index[i] = row_index[i];
    }
    time = omp_get_wtime() - start;
    return 0;
}

// Task 2
double CPUSpeed()
{
    wchar_t Buffer[_MAX_PATH];
    DWORD BufSize = _MAX_PATH;
    DWORD dwMHz = _MAX_PATH;
    HKEY hKey;

    // open the key where the proc speed is hidden:
    long lError = RegOpenKeyExW(HKEY_LOCAL_MACHINE,
                                L"HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",
                                0,
                                KEY_READ,
                                &hKey);
    if(lError != ERROR_SUCCESS)
    {// if the key is not found, tell the user why:
        FormatMessageW(FORMAT_MESSAGE_FROM_SYSTEM,
                       NULL,
                       lError,
                       0,
                       Buffer,
                       _MAX_PATH,
                       0);
        wprintf(Buffer);
        return 0;
    }

    // query the key:
    RegQueryValueExW(hKey, L"~MHz", NULL, NULL, (LPBYTE) &dwMHz, &BufSize);
    return (double)dwMHz;
}

void get_cpu_info(){
    char CPUBrandString[0x40];
    unsigned int CPUInfo[4] = {0,0,0,0};

    __cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    unsigned int nExIds = CPUInfo[0];

    memset(CPUBrandString, 0, sizeof(CPUBrandString));

    for (unsigned int i = 0x80000000; i <= nExIds; ++i)
    {
        __cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);

        if (i == 0x80000002)
            memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000003)
            memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000004)
            memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
    }

    std::cout << "CPU Type: " << CPUBrandString << std::endl;
    std::cout << "Number of cores: " << omp_get_num_procs() << std::endl;
    std::cout << "CPU Speed: " << CPUSpeed() / 1000 << "GHz" << std::endl;
    MEMORYSTATUSEX statex;
    statex.dwLength = sizeof (statex);
    GlobalMemoryStatusEx (&statex);
    std::cout << "Memory: " << (float)statex.ullTotalPhys/(1024*1024) << "MB" << std::endl;
}

void get_GPU_info(){
    int i, j;
    char* value;
    size_t valueSize;
    cl_uint platformCount;
    cl_platform_id* platforms;
    cl_uint deviceCount;
    cl_device_id* devices;
    cl_uint maxComputeUnits;

    // get all platforms
    clGetPlatformIDs(0, NULL, &platformCount);
    platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id) * platformCount);
    clGetPlatformIDs(platformCount, platforms, NULL);

    for (i = 0; i < platformCount; i++) {

        // get all devices
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
        devices = (cl_device_id*)malloc(sizeof(cl_device_id) * deviceCount);
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);

        // for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {

            // print device name
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
            value = (char*)malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
            printf("%d. Device: %s\n", j + 1, value);
            free(value);

            // print hardware device version
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
            value = (char*)malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
            printf(" %d.%d Hardware version: %s\n", j + 1, 1, value);
            free(value);

            // print software driver version
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
            value = (char*)malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
            printf(" %d.%d Software version: %s\n", j + 1, 2, value);
            free(value);

            // print c version supported by compiler for device
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
            value = (char*)malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            printf(" %d.%d OpenCL C version: %s\n", j + 1, 3, value);
            free(value);

            // print parallel compute units
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                            sizeof(maxComputeUnits), &maxComputeUnits, NULL);
            printf(" %d.%d Parallel compute units: %d\n", j + 1, 4, maxComputeUnits);

        }

        free(devices);

    }

    free(platforms);
}

// Show gnu graphs
template<size_t cols>
void show_graph(int task_num, double* x, double (&y)[3][cols], int size, std::string title){
    FILE* fp;
    std::string file_name = "task_" + std::to_string(task_num) + "_data.tmp";
    fp = std::fopen(file_name.c_str(), "w");
    FILE* gnupipe;
    gnupipe = _popen("gnuplot -persistent", "w");
    for(int i = 0; i < size; i++){
        fprintf(fp, "%f %f %f %f\n", x[i], y[0][i], y[1][i], y[2][i]);
    }


    std::fprintf(gnupipe, ("set title \"" + title + "\"\n").c_str());

    fprintf(gnupipe, "set xlabel \"Matrix size\"\n"
                     "set ylabel \"Time (s)\"\n");

    std::fprintf(gnupipe, ("plot \"" + file_name + "\" u 1:2 with lines title \"Threads = 2\",\\\n"
                     "\"" + file_name + "\" u 1:3 with lines title \"Threads = 4\",\\\n"
                     "\"" + file_name + "\" u 1:4 with lines title \"Threads = 8\",\\\n").c_str());
    fprintf(gnupipe, "exit");
}


int main() {
    // Task 1
    printf("Dinislam Gabitov. Programming assignment #1\n");
    get_cpu_info();
    //get_GPU_info();

    /*
    * ------------------------------- Task 11
    */

    crsMatrix task_11_crs;
    regMatrix task_11_reg;
    double *task_11_v1, *task_11_out1, *task_11_out2;
    double task_11_diff;
    double task_11_time1, task_11_time2;
    std::cout << "Task 11 - CRS-vector mult check" << std::endl;
    for(int i = 10; i <= 1000; i = i * 10){
        AllocateCRSMatrix(i, i * 3, task_11_crs);
        GenerateSpecialCRS(1, i, 2, task_11_crs);
        CRSToRegular(task_11_crs, task_11_reg);
        AllocateVector(i, &task_11_v1);
        AllocateVector(i, &task_11_out1);
        AllocateVector(i, &task_11_out2);

        for(int j = 0; j < i; j++){
            task_11_v1[j] = j;
        }
        SeqCRSMult(task_11_crs, task_11_v1, task_11_out1, task_11_time1);
        SeqRegMult(task_11_reg, task_11_v1, task_11_out2, task_11_time2);
        CompareVectors(task_11_out1, task_11_out2, i, task_11_diff);
        printf("For n = %d ", i);
        if(task_11_diff < 0.05){
            printf("Correct\n");
        } else{
            printf("Incorrect\n");
        }
        FreeCRSMatrix(task_11_crs);
        FreeRegMatrix(task_11_reg);
        FreeVector(&task_11_v1);
        FreeVector(&task_11_out1);
        FreeVector(&task_11_out2);
    }

    /*
     * ------------------------------- Task 13
     */

    crsMatrix task_13_crs;
    regMatrix task_13_reg;
    double *task_13_v1, *task_13_out;
    double task_13_time;
    double task_13_total_time[3][10];
    double task_13_x_plot[10];
    for(int i = 500; i <= 5000; i += 500){
        task_13_x_plot[i / 500 - 1] = i;

        AllocateCRSMatrix(i, i * 3, task_13_crs);
        GenerateSpecialCRS(1, i, i / 100, task_13_crs);
        CRSToRegular(task_13_crs, task_13_reg);
        AllocateVector(i, &task_13_v1);
        AllocateVector(i, &task_13_out);
        for(int j = 0; j < i; j++){
            task_13_v1[j] = j;
        }

        OmpRegMult(task_13_reg, task_13_v1, task_13_out, 2, task_13_time);
        task_13_total_time[0][i / 500 - 1] = task_13_time;

        OmpRegMult(task_13_reg, task_13_v1, task_13_out, 4, task_13_time);
        task_13_total_time[1][i / 500 - 1] = task_13_time;

        OmpRegMult(task_13_reg, task_13_v1, task_13_out, 8, task_13_time);
        task_13_total_time[2][i / 500 - 1] = task_13_time;

        FreeCRSMatrix(task_13_crs);
        FreeRegMatrix(task_13_reg);
        FreeVector(&task_13_v1);
        FreeVector(&task_13_out);
    }

    show_graph(13, task_13_x_plot, task_13_total_time, 10, "Regular matrix-vector mult");

    /*
     * ------------------------------- Task 15
     */

    crsMatrix task_15_crs;
    double *task_15_v1, *task_15_out;
    double task_15_time;
    double task_15_total_time[3][8];
    double task_15_x_plot[8];
    for(int i = 5000; i <= 40000; i += 5000){
        task_16_x_plot[i / 5000 - 1] = i;

        AllocateCRSMatrix(i, i * 3, task_15_crs);
        GenerateSpecialCRS(1, i, i / 100, task_15_crs);
        AllocateVector(i, &task_15_v1);
        AllocateVector(i, &task_15_out);
        for(int j = 0; j < i; j++){
            task_15_v1[j] = j;
        }

        OmpCRSMult(task_15_crs, task_15_v1, task_15_out, 2, task_15_time);
        task_15_total_time[0][i / 5000 - 1] = task_15_time;

        OmpCRSMult(task_15_crs, task_15_v1, task_15_out, 4, task_15_time);
        task_15_total_time[1][i / 5000 - 1] = task_15_time;

        OmpCRSMult(task_15_crs, task_15_v1, task_15_out, 8, task_15_time);
        task_15_total_time[2][i / 5000 - 1] = task_15_time;

        FreeCRSMatrix(task_15_crs);
        FreeVector(&task_15_v1);
        FreeVector(&task_15_out);
    }

    show_graph(15, task_15_x_plot, task_15_total_time, 8, "CRS matrix-vector mult");

    /*
     * ---------------------------- Task 18
     */
    std::cout << "Task 18 - CRS matrix-matrix mult check" << std::endl;

    crsMatrix task_18_crs_1, task_18_crs_2, task_18_crs_out;
    regMatrix task_18_reg_1, task_18_reg_2, task_18_reg_out, task_18_crs_to_reg;

    double task_18_diff, task_18_time;
    for(int i = 10; i <= 1000; i *= 10){
        AllocateCRSMatrix(i, i * 3, task_18_crs_1);
        GenerateSpecialCRS(1, i, 2, task_18_crs_1);
        CRSToRegular(task_18_crs_1, task_18_reg_1);

        AllocateCRSMatrix(i, i * 3, task_18_crs_2);
        GenerateSpecialCRS(1, i, 2, task_18_crs_2);
        CRSToRegular(task_18_crs_2, task_18_reg_2);

        SeqCRSMatMult(task_18_crs_1, task_18_crs_2, task_18_crs_out, task_18_time);
        SeqRegMatMult(task_18_reg_1, task_18_reg_2, task_18_reg_out, task_18_time);

        CRSToRegular(task_18_crs_out, task_18_crs_to_reg);
        CompareRegMatrix(task_18_reg_out, task_18_crs_to_reg, task_18_diff);

        printf("For n = %d ", i);
        if(task_18_diff < 0.005){
            printf("Correct\n");
        } else{
            printf("Incorrect\n");
            std::cout << task_18_diff << std::endl;
        }
        FreeCRSMatrix(task_18_crs_1);
        FreeCRSMatrix(task_18_crs_2);
        FreeCRSMatrix(task_18_crs_out);
        FreeRegMatrix(task_18_reg_1);
        FreeRegMatrix(task_18_reg_2);
        FreeRegMatrix(task_18_reg_out);
        FreeRegMatrix(task_18_crs_to_reg);

    }

    /*
     * ---------------------------- Task 20
     */

    crsMatrix task_20_crs_1, task_20_crs_2;
    regMatrix task_20_reg_1, task_20_reg_2, task_20_reg_out;
    double task_20_time;
    double task_20_total_time[3][20];
    double task_20_x_plot[20];

    for(int i = 50; i <= 1000; i += 50){
        task_20_x_plot[i / 50 - 1] = i;

        AllocateCRSMatrix(i, i * 3, task_20_crs_1);
        GenerateSpecialCRS(i, i, i / 10, task_20_crs_1);
        AllocateCRSMatrix(i, i * 3, task_20_crs_2);
        GenerateSpecialCRS(i + 1, i, i / 10, task_20_crs_2);

        AllocateRegMatrix(i, task_20_reg_1);
        CRSToRegular(task_20_crs_1, task_20_reg_1);
        AllocateRegMatrix(i, task_20_reg_2);
        CRSToRegular(task_20_crs_2, task_20_reg_2);
        AllocateRegMatrix(i, task_20_reg_out);


        OmpRegMatMult(task_20_reg_1, task_20_reg_2, task_20_reg_out, 2, task_20_time);
        task_20_total_time[0][i / 50 - 1] = task_20_time;

        OmpRegMatMult(task_20_reg_1, task_20_reg_2, task_20_reg_out, 4, task_20_time);
        task_20_total_time[1][i / 50 - 1] = task_20_time;

        OmpRegMatMult(task_20_reg_1, task_20_reg_2, task_20_reg_out, 8, task_20_time);
        task_20_total_time[2][i / 50 - 1] = task_20_time;

        FreeCRSMatrix(task_20_crs_1);
        FreeCRSMatrix(task_20_crs_2);
        FreeRegMatrix(task_20_reg_1);
        FreeRegMatrix(task_20_reg_2);
        FreeRegMatrix(task_20_reg_out);
    }


    show_graph(20, task_20_x_plot, task_20_total_time, 20, "Regular matrix-matrix mult");


    /*
     * ---------------------------- Task 22
     */

    crsMatrix task_22_crs_1, task_22_crs_2, task_22_crs_out;
    double task_22_time;
    double task_22_total_time[3][20];
    double task_22_x_plot[20];

    for(int i = 50; i <= 1000; i += 50){
        task_22_x_plot[i / 50 - 1] = i;

        AllocateCRSMatrix(i, i * 3, task_22_crs_1);
        GenerateSpecialCRS(i, i, i / 20, task_22_crs_1);
        AllocateCRSMatrix(i, i * 3, task_22_crs_2);
        GenerateSpecialCRS(i + 1, i, i / 20, task_22_crs_2);

        OmpCRSMatMult(task_22_crs_1, task_22_crs_2, task_22_crs_out, 2, task_22_time);
        task_22_total_time[0][i / 50 - 1] = task_22_time;
        FreeCRSMatrix(task_22_crs_out);

        OmpCRSMatMult(task_22_crs_1, task_22_crs_2, task_22_crs_out, 4, task_22_time);
        task_22_total_time[1][i / 50 - 1] = task_22_time;
        FreeCRSMatrix(task_22_crs_out);

        OmpCRSMatMult(task_22_crs_1, task_22_crs_2, task_22_crs_out, 8, task_22_time);
        task_22_total_time[2][i / 50 - 1] = task_22_time;
        FreeCRSMatrix(task_22_crs_out);

        FreeCRSMatrix(task_22_crs_1);
        FreeCRSMatrix(task_22_crs_2);
    }

    show_graph(22, task_22_x_plot, task_22_total_time, 20, "CRS matrix-matrix mult");

    return 0;

}