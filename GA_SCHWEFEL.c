/*GA_SCHWEFEL*/
/*
程式會輸出兩個檔案
1."GA_result.txt", 記錄三種指定的fitness function平均值
2."GA_chromosome.txt", 記錄第一次repeat時的所有點座標(**注意此檔案是用'寫入在檔案後'方式，所以每次執行前記得先
刪除該檔案，避免重複寫入同檔案)

另由matlab繪製圖表
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define POPULATION_SIZE 96		//指定 96
#define CHROMOSOME_LENGTH 24	//指定 24
#define CROSSOVER_RATE 0.7
#define MUTATE_RATE 0.01
#define ITERATION 20			//自訂20
#define REPEAT 10				
#define MAXIMUM 500				//測試範圍
#define MINIMUM -500
#define SCORE_TYPE 3			// 0:total average, 1:top-48 average, 2:top-20 average

void initialGroup(int arr[POPULATION_SIZE][CHROMOSOME_LENGTH]) {
	int i, j;
	int num;
	int* random = &num;

	for (i = 0; i < POPULATION_SIZE; i++) {
		for (j = 0; j < CHROMOSOME_LENGTH; j++) {
			*random = rand();
			(*random > RAND_MAX / 2) ? (arr[i][j] = 1) : (arr[i][j] = 0);
		}
	}
}
double binToDec(int arr[]) {
	double num = 0;
	int i;
	for (i = 0; i < CHROMOSOME_LENGTH / 2; i++) {
		//printf("%d", arr[i]);
		num += arr[(CHROMOSOME_LENGTH / 2 -1)-i] * pow(2, i);
	}
	//printf("toDec:%lf\n", num);
	return num;
}
void fitnessSHCWEFEL(int inputArr[POPULATION_SIZE][CHROMOSOME_LENGTH], double outputFitnessArr[POPULATION_SIZE]) {
	int x1_bin[CHROMOSOME_LENGTH/2] = { 0 }, x2_bin[CHROMOSOME_LENGTH/2] = { 0 };
	double x1_dec, x2_dec;
	int i,j;

	for (i = 0; i < POPULATION_SIZE; i++) {
		//BIN to DEC
		for (j = 0; j < CHROMOSOME_LENGTH / 2; j++) {
			x1_bin[j] = inputArr[i][j];
		}
		for (j = CHROMOSOME_LENGTH / 2; j < CHROMOSOME_LENGTH; j++) {
			x2_bin[j - CHROMOSOME_LENGTH / 2] = inputArr[i][j];
		}
		x1_dec = binToDec(x1_bin);
		x2_dec = binToDec(x2_bin);
		//normalized to [0:1]
		x1_dec = (x1_dec - 0) / ((pow(2, CHROMOSOME_LENGTH / 2) - 1) - 0);
		x2_dec = (x2_dec - 0) / ((pow(2, CHROMOSOME_LENGTH / 2) - 1) - 0);
		//transfer to our scope
		x1_dec = x1_dec * (MAXIMUM - MINIMUM) + MINIMUM; 
		x2_dec = x2_dec * (MAXIMUM - MINIMUM) + MINIMUM;
		//calculate fitness score
		if (x1_dec < 0) x1_dec = -x1_dec;		//abs()
		if (x2_dec < 0) x2_dec = -x2_dec;
		outputFitnessArr[i] = 1 / (418.9829 * 2 - (x1_dec * sin(sqrt(x1_dec)) + x2_dec * sin(sqrt(x2_dec))));
		////1 / (418.9829 * 2 - ( x1_dec * sin(sqrt(x1_dec)) + x2_dec * sin(sqrt(x2_dec)) ) )

		////printf("(%4.2lf, %4.2lf) = %lf\n", x1_dec, x2_dec, outputFitnessArr[i]);
	}
}	
void sortDoubleArr(double arr[POPULATION_SIZE], int index[POPULATION_SIZE]) {
	int i,j,indexTemp, Flag=1;
	double temp;

	for (i = 0; i < POPULATION_SIZE; i++) {						//fill index
		index[i] = i;
	}
	for (i = 0; i < POPULATION_SIZE-1 && Flag; i++) {			//sorting
		Flag = 0;
		for (j = 0; j < POPULATION_SIZE-1-i; j++) {
			if (arr[j] > arr[j+1]) {
				temp = arr[j];
				arr[j] = arr[j+1];
				arr[j+1] = temp;

				indexTemp = index[j];							//紀錄排序前的位置
				index[j] = index[j + 1];
				index[j + 1] = indexTemp;

				Flag = 1;
			}
		}
	}
}
void selectAndCopy(int parentArr[POPULATION_SIZE][CHROMOSOME_LENGTH], int childArr[POPULATION_SIZE][CHROMOSOME_LENGTH], double fitnessArr[POPULATION_SIZE]) {
	int i, j, k, bestIndex, currentIndex, sortIndex[POPULATION_SIZE] = { 0 };
	double selectIndex, fitnessSum = 0, currentSum;

	for (i = 0; i < POPULATION_SIZE; i++) {
		fitnessSum += fitnessArr[i];
	}
	/*
		for (i = 0; i < POPULATION_SIZE; i++) {
		printf("%lf ", fitnessArr[i]);
		printf("\n");
	}
	*/
	sortDoubleArr(fitnessArr, sortIndex);
	/*
		for (i = 0; i < POPULATION_SIZE; i++) {
		printf("%d ", sortIndex[i]+1);
		printf("\n");
	}
	for (i = 0; i < POPULATION_SIZE; i++) {
		printf("%lf ", fitnessArr[i]);
		printf("\n");
	}
	*/

	//save elite
	bestIndex = sortIndex[POPULATION_SIZE - 1];
	printf("best score(取倒數):%lf\n", fitnessArr[POPULATION_SIZE - 1]);
	for (i = 0; i < CHROMOSOME_LENGTH; i++) {
		childArr[0][i] = parentArr[bestIndex][i];
	}

	//random select and copy
	for (i = 1; i < POPULATION_SIZE; i++) {
		selectIndex = (double)rand()/RAND_MAX * fitnessSum;
		currentSum = 0;
		j = 0;
		////printf("i=%d, selectIndex=%lf\n", i, selectIndex);
		while (j < POPULATION_SIZE) {							//累計輪盤
			currentSum += fitnessArr[j];
			if (selectIndex <= currentSum) {
				currentIndex = sortIndex[j];					//取出排序前的父代位置
				//printf("父%d => 子%d\n", j+1, i+1);
				for (k = 0; k < CHROMOSOME_LENGTH; k++) {
					childArr[i][k] = parentArr[currentIndex][k];
				}
				break;
			}
			j++;
		}
	}
}
void crossover_1point(int arr[POPULATION_SIZE][CHROMOSOME_LENGTH], int rowIndex){
	int i;
	int crossIndex = rand()%CHROMOSOME_LENGTH;
	//printf("crossIndex=%d\n", crossIndex+1);
	for (i = crossIndex; i < CHROMOSOME_LENGTH; i++) {
		(arr[rowIndex][i] == 0) ? (arr[rowIndex][i] = 1) : (arr[rowIndex][i] = 0);
	}
}
void crossover_2point(int arr[POPULATION_SIZE][CHROMOSOME_LENGTH], int rowIndex) {
	int i,temp;
	int startCrossIndex = rand() % CHROMOSOME_LENGTH;
	int endCrossIndex = rand() % CHROMOSOME_LENGTH;

	if (startCrossIndex > endCrossIndex) {
		temp = startCrossIndex;
		startCrossIndex = endCrossIndex;
		endCrossIndex = temp;
		//printf("exchange!!\n");
	}
	//printf("startCrossIndex=%d, endCrossIndex=%d\n", startCrossIndex+1,endCrossIndex + 1);
	for (i = startCrossIndex; i <= endCrossIndex; i++) {
		(arr[rowIndex][i] == 0) ? (arr[rowIndex][i] = 1) : (arr[rowIndex][i] = 0);
	}
}
void mutate_1point(int arr[POPULATION_SIZE][CHROMOSOME_LENGTH], int rowIndex) {
	int index = rand() % CHROMOSOME_LENGTH;
	
	//printf("mutate:%d\n", index+1);
	(arr[rowIndex][index] == 0) ? (arr[rowIndex][index] = 1) : (arr[rowIndex][index] = 0);
}
void recordFitness(int inputArr[POPULATION_SIZE][CHROMOSOME_LENGTH], double resultArr[SCORE_TYPE][REPEAT][ITERATION], int repeat, int iteration) {
	int i, j;
	int x1_bin[CHROMOSOME_LENGTH / 2] = { 0 }, x2_bin[CHROMOSOME_LENGTH / 2] = { 0 };
	int sortIndex[POPULATION_SIZE] = { 0 };
	double x1_dec, x2_dec, fitnessSum = 0;
	double fitnessArr[POPULATION_SIZE] = { 0 };

	/* 算出原版的fitness(不取倒數) */
	for (i = 0; i < POPULATION_SIZE; i++) {
		//BIN to DEC
		for (j = 0; j < CHROMOSOME_LENGTH / 2; j++) {
			x1_bin[j] = inputArr[i][j];
		}
		for (j = CHROMOSOME_LENGTH / 2; j < CHROMOSOME_LENGTH; j++) {
			x2_bin[j - CHROMOSOME_LENGTH / 2] = inputArr[i][j];
		}
		x1_dec = binToDec(x1_bin);
		x2_dec = binToDec(x2_bin);
		//normalized to [0:1]
		x1_dec = (x1_dec - 0) / ((pow(2, CHROMOSOME_LENGTH / 2) - 1) - 0);
		x2_dec = (x2_dec - 0) / ((pow(2, CHROMOSOME_LENGTH / 2) - 1) - 0);
		//transfer to our scope
		x1_dec = x1_dec * (MAXIMUM - MINIMUM) + MINIMUM;
		x2_dec = x2_dec * (MAXIMUM - MINIMUM) + MINIMUM;
		//calculate fitness score
		if (x1_dec < 0) x1_dec = -x1_dec;		//abs()
		if (x2_dec < 0) x2_dec = -x2_dec;
		fitnessArr[i] = 418.9829 * 2 - (x1_dec * sin(sqrt(x1_dec)) + x2_dec * sin(sqrt(x2_dec)));
	}

	for (i = 0; i < POPULATION_SIZE; i++) {
		fitnessSum += fitnessArr[i];
	}
	sortDoubleArr(fitnessArr, sortIndex);

	//total average
	resultArr[0][repeat][iteration] = fitnessSum / POPULATION_SIZE;

	//top-48 average
	fitnessSum = 0;
	for (i = 0; i < 48; i++) {
		fitnessSum += fitnessArr[i];
	}
	resultArr[1][repeat][iteration] = fitnessSum / 48;

	//top-20 average
	fitnessSum = 0;
	for (i = 0; i < 20; i++) {
		fitnessSum += fitnessArr[i];
	}
	resultArr[2][repeat][iteration] = fitnessSum / 20;
}
void randomSortAndCopy(int inputArr[POPULATION_SIZE][CHROMOSOME_LENGTH], int outputArr[POPULATION_SIZE][CHROMOSOME_LENGTH]) {
	int i, j, index, tempArr[CHROMOSOME_LENGTH] = { 0 };

	for (i = 0; i < POPULATION_SIZE; i++) {				//copy
		for (j = 0; j < CHROMOSOME_LENGTH; j++) {
			outputArr[i][j] = inputArr[i][j];
		}
	}
	for (i = 0; i < POPULATION_SIZE; i++) {				//random sort
		index = rand() % POPULATION_SIZE;
		//printf("%d換到%d\n", index + 1, i+1);
		for (j = 0; j < CHROMOSOME_LENGTH; j++) tempArr[j] = outputArr[i][j];
		for (j = 0; j < CHROMOSOME_LENGTH; j++) outputArr[i][j] = outputArr[index][j];
		for (j = 0; j < CHROMOSOME_LENGTH; j++) outputArr[index][j] = tempArr[j];
	}
}
void outputFile(double arr[SCORE_TYPE][REPEAT][ITERATION]) {
	int i, j, k;
	FILE* fptr;						

	fopen_s(&fptr, "GA_result.txt", "w");

	if (fptr != NULL) {		
		for (i = 0; i < SCORE_TYPE; i++) {
			switch (i) {
			case 0: fprintf(fptr, "total average:\n"); break;
			case 1: fprintf(fptr, "top-48 average:\n"); break;
			case 2: fprintf(fptr, "top-20 average:\n"); break;
			}
			for (j = 0; j < REPEAT; j++) {
				for (k = 0; k < ITERATION; k++) {
					fprintf(fptr, "%lf ", arr[i][j][k]);
				}
				fprintf(fptr,"\n");
			}
		}
		fclose(fptr);
		printf("檔案輸出完成\n");
	}
	else
		printf("檔案開啟失敗\n");
}
void outputChromosome(int inputArr[POPULATION_SIZE][CHROMOSOME_LENGTH],int iteration) {
	int i, j;
	int x1_bin[CHROMOSOME_LENGTH / 2] = { 0 }, x2_bin[CHROMOSOME_LENGTH / 2] = { 0 };
	double x1_dec, x2_dec;
	FILE* fptr;

	fopen_s(&fptr, "GA_chromosome.txt", "a");

	if (fptr != NULL) {
		fprintf(fptr, "----------ITERATION:%d----------\n",iteration+1);

		for (i = 0; i < POPULATION_SIZE; i++) {
			//BIN to DEC
			for (j = 0; j < CHROMOSOME_LENGTH / 2; j++) {
				x1_bin[j] = inputArr[i][j];
			}
			for (j = CHROMOSOME_LENGTH / 2; j < CHROMOSOME_LENGTH; j++) {
				x2_bin[j - CHROMOSOME_LENGTH / 2] = inputArr[i][j];
			}
			x1_dec = binToDec(x1_bin);
			x2_dec = binToDec(x2_bin);
			//normalized to [0:1]
			x1_dec = (x1_dec - 0) / ((pow(2, CHROMOSOME_LENGTH / 2) - 1) - 0);
			x2_dec = (x2_dec - 0) / ((pow(2, CHROMOSOME_LENGTH / 2) - 1) - 0);
			//transfer to our scope
			x1_dec = x1_dec * (MAXIMUM - MINIMUM) + MINIMUM;
			x2_dec = x2_dec * (MAXIMUM - MINIMUM) + MINIMUM;

			fprintf(fptr, "(%4.2lf,%4.2lf) ", x1_dec, x2_dec);
		}
		fprintf(fptr, "\n");
		fclose(fptr);
	//	printf("檔案輸出完成\n");
	}
	//else
	//	printf("檔案開啟失敗\n");
}

void printChromosome(int arr[POPULATION_SIZE][CHROMOSOME_LENGTH]) {
	int i, j;
	for (i = 0; i < POPULATION_SIZE; i++) {
		for (j = 0; j < CHROMOSOME_LENGTH; j++) {
			printf("%d ", arr[i][j]);
		}
		printf("\n");
	}
}
void printFitness(double arr[POPULATION_SIZE]) {
	int i;
	for (i = 0; i < POPULATION_SIZE; i++) {
			printf("%lf ", arr[i]);
		printf("\n");
	}
}
int main(void) {
	int chromosome_parent[POPULATION_SIZE][CHROMOSOME_LENGTH] = { 0 };
	int chromosome_child[POPULATION_SIZE][CHROMOSOME_LENGTH] = { 0 };
	double fitnessArray[POPULATION_SIZE] = { 0 };
	double result[SCORE_TYPE][REPEAT][ITERATION] = { 0 };

	int repeat, iteration, i;
	srand(time(NULL));

	//Genitic Algorithm
	for (repeat = 0; repeat < REPEAT; repeat++) {					//重複次數
		initialGroup(chromosome_parent);
		for (iteration = 0; iteration < ITERATION; iteration++) {	//迭代次數
			fitnessSHCWEFEL(chromosome_parent, fitnessArray);

			printf("repeat:%d, iteration:%d\n", repeat+1, iteration+1);
			////printChromosome(chromosome_child);

			//select and copy
			selectAndCopy(chromosome_parent, chromosome_child, fitnessArray);

			//crossover
			for (i = 0; i < POPULATION_SIZE; i++) {
				if ((double)rand() / RAND_MAX <= CROSSOVER_RATE) {
					crossover_2point(chromosome_child, i);
					////printf( "%d CROSSOVER %lf\n", i, (double)rand() / RAND_MAX);
				}
			}

			//mutate
			for (i = 0; i < POPULATION_SIZE; i++) {
				if ((double)rand() / RAND_MAX <= MUTATE_RATE) {
					mutate_1point(chromosome_child, i);
					////printf("%d MUTATE %lf\n", i, (double)rand() / RAND_MAX);
				}
			}

			//record and reload
			fitnessSHCWEFEL(chromosome_child, fitnessArray);
			recordFitness(chromosome_child, result, repeat, iteration);

			if(repeat==0) outputChromosome(chromosome_child, iteration);

			if (result[0][repeat][iteration] < 1e-5) break;		//如果已經收斂，結束此次迭代
			randomSortAndCopy(chromosome_child, chromosome_parent);
			
			printf("--------------------------\n");
		}
		
	}
	outputFile(result);

	/*
		printf("%d after\n", i+1);
		printChromosome(chromosome_child);

	*/
	return 0;
}