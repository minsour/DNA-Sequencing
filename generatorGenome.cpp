#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <random>
#include <fstream>
#include <list>
#include <map>
#include <chrono>
using namespace std;
#define FILE_SIZE 1
#define N 1000000
#define M 50000
#define TOTAL_LEN N*FILE_SIZE
#define DIFF_PERCENT 0.01
#define L 32
#define ALLOW_MISMATCHED 2

const char DNA_LIST[] = { "ACGT" };
const string folderName[] = { "reference500000", "my500000", "short500000", "recovery500000" };
int changeIndex[int(TOTAL_LEN * DIFF_PERCENT)];
set<unsigned long long int> shortRead;
unsigned long long int bitGenome[N];
vector<int> mismatched[FILE_SIZE];
bool matched[FILE_SIZE][N];
int shortReadIndex[M]; 

string createFilePath(int fileType, int fileIndex) {
	return folderName[fileType] + "\\genome" + to_string(fileIndex) + ".txt";
}

int getRandomNumber(int min, int max) {
	//< 1�ܰ�. �õ� ����
	random_device rn;
	mt19937_64 rnd(rn());

	//< 2�ܰ�. ���� ���� ( ���� )
	uniform_int_distribution<int> range(min, max);

	//< 3�ܰ�. �� ����
	return range(rnd);
}

unsigned char toBits(char dna) {
	switch (dna) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		break;
	}
}

void makeReferenceGenome() {
	for (int i = 0; i < FILE_SIZE; i++) {
		string referenceGenomeSeq;
		cout << "���� " << N << "�� ReferenceGenomeSeq" << i << " ���� ��..!" << endl;
		for (int j = 0; j < N; ++j) {
			referenceGenomeSeq += DNA_LIST[getRandomNumber(0, 3)];
		}
		ofstream *referenceGenomeFile = new ofstream(createFilePath(0, i));
		if (referenceGenomeFile->is_open()) {
			*referenceGenomeFile << referenceGenomeSeq;
			referenceGenomeFile->close();
			delete referenceGenomeFile;
		}
		else {
			cout << "ReferenceGenome File ���� �Ұ�" << endl;
			exit(1);
		}
		cout << "ReferenceGenomeSeq" << i << " ���� �Ϸ�!" << endl;
	}
}

void makeMyGenome() {
	cout << "MyGenome ���� Index ���� ��..!" << endl;
	map <int, bool> checkDuplication;
	int diffSize = TOTAL_LEN * DIFF_PERCENT;
	for (int i = 0; i < diffSize; i++) {
		int randomNumber;
		while (true) {
			randomNumber = getRandomNumber(0, TOTAL_LEN);
			if (checkDuplication.count(randomNumber) <= 0) {
				break;
			}
		}
		checkDuplication[randomNumber] = true;
		changeIndex[i] = randomNumber;
	}

	sort(changeIndex, changeIndex + diffSize);
	cout << "MyGenome ���� Index ���� �Ϸ�..!" << endl;

	int index = 0;
	for (int i = 0; i < FILE_SIZE && index < diffSize; i++) {
		cout << i << "��° ReferenceGenome�� MyGenome ���� ��..!" << endl;

		ifstream *referenceGenomeFile = new ifstream(createFilePath(0, i));
		string genome;
		if (referenceGenomeFile->is_open()) {
			*referenceGenomeFile >> genome;
			referenceGenomeFile->close();
			delete referenceGenomeFile;
		}
		else {
			cout << "ReferenceGenome File �б� �Ұ�" << endl;
			exit(1);
		}

		if (changeIndex[index] < i * N || changeIndex[index] >= (i + 1) * N) {
			ofstream *myGenomeFile = new ofstream(createFilePath(1, i));
			if (myGenomeFile->is_open()) {
				*myGenomeFile << genome;
				myGenomeFile->close();
				delete myGenomeFile;
			}
			else {
				cout << "MyGenome File ���� �Ұ�" << endl;
				exit(1);
			}
			cout << i << "��° ReferenceGenome�� MyGenome ���� �Ϸ�..!" << endl;
			continue;
		}

		while (index < diffSize && changeIndex[index] >= i * N && changeIndex[index] < (i + 1) * N) {
			char temp = genome[changeIndex[index] - (i * N)];
			char changeGenome;
			while ((changeGenome = DNA_LIST[getRandomNumber(0, 3)]) == temp);
			genome[changeIndex[index] - (i * N)] = changeGenome;
			index++;
		}

		ofstream *myGenomeFile = new ofstream(createFilePath(1, i));
		if (myGenomeFile->is_open()) {
			*myGenomeFile << genome;
			myGenomeFile->close();
			delete myGenomeFile;
		}
		else {
			cout << "MyGenome File ���� �Ұ�" << endl;
			exit(1);
		}

		cout << i << "��° ReferenceGenome�� MyGenome ���� �Ϸ�..!" << endl;
	}
}

void makeShortReadGenome() {
	cout << "ShortReadGenome Index ���� ��..!" << endl;

	for (int i = 0; i < M; i++) {
		shortReadIndex[i] = getRandomNumber(0, TOTAL_LEN - L);
	}

	sort(shortReadIndex, shortReadIndex + M);
	cout << "ShortReadGenome Index ���� �Ϸ�..!" << endl;

	int index = 0;
	ofstream *shortReadFile = new ofstream(createFilePath(2, 0));
	if (shortReadFile->is_open() == false) {
		cout << "Short Read File ���� �Ұ�" << endl;
		exit(1);
	}

	for (int i = 0; i < FILE_SIZE && index < M; i++) {
		if (shortReadIndex[index] < i * N || shortReadIndex[index] >= (i + 1) * N) {
			continue;
		}
		ifstream *referenceGenomeFile = new ifstream(createFilePath(1, i));
		string genome;
		if (referenceGenomeFile->is_open()) {
			*referenceGenomeFile >> genome;
			referenceGenomeFile->close();
			delete referenceGenomeFile;
		}
		else {
			cout << "ReferenceGenome File �б� �Ұ�" << endl;
			exit(1);
		}

		while (index < M && shortReadIndex[index] >= i * N && shortReadIndex[index] < (i + 1) * N) {
			if (shortReadIndex[index] - (i * N) + L > genome.length()) {
				referenceGenomeFile = new ifstream(createFilePath(1, i + 1));
				string temp;
				if (referenceGenomeFile->is_open()) {
					*referenceGenomeFile >> temp;
					referenceGenomeFile->close();
					delete referenceGenomeFile;
				}
				else {
					cout << "ReferenceGenome File �б� �Ұ�" << endl;
					exit(1);
				}
				genome += temp;
			}
			*shortReadFile << genome.substr(shortReadIndex[index] - i * N, L) << endl;
			index++;
		}
	}
}

int main() {
	makeReferenceGenome();
	makeMyGenome();
	makeShortReadGenome();
	return 0;
}