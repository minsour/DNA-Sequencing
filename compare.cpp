/* 2013113506 �ڹμ� */
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <random>
#include <fstream>
#include <list>
#include <map>
#include <chrono>
using namespace std;

#define FILE_SIZE 1
#define N 500000
#define M 50000
#define TOTAL_LEN N*FILE_SIZE
#define DIFF_PERCENT 0.01
#define L 32
#define ALLOW_MISMATCHED 2

const char DNA_LIST[] = { "ACGT" };
const string folderName[] = { "reference500000", "my500000", "short500000", "recovery" };
int changeIndex[int(TOTAL_LEN * DIFF_PERCENT)];
multiset<unsigned long long int> shortRead;
unsigned long long int bitGenome[N];
vector<int> mismatched[FILE_SIZE];
bool matched[FILE_SIZE][N];

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

void readShortRead() {
	ifstream *shortReadGenomeFile = new ifstream(createFilePath(2, 0));
	if (shortReadGenomeFile->is_open()) {
		for (int i = 0; i < M; i++) {
			if (i % 1000000 == 0) {
				cout << i << "/" << M << " ShortRead �д� ��..." << endl;
			}
			string temp;
			*shortReadGenomeFile >> temp;
			unsigned long long int data = 0x0;
			for (int j = 0; j < L; j++) {
				data <<= 2;
				data |= toBits(temp[j]);
			}
			shortRead.insert(data);
		}
		shortReadGenomeFile->close();
		delete shortReadGenomeFile;
	}
	else {
		cout << "ShortRead Genome File �б� �Ұ�" << endl;
		exit(1);
	}
}

void recoveryMatchedGenome() {
	string fixedMyGenome;
	string nextReferenceGenomeSeq;
	ifstream *referenceGenomeFile = new ifstream(createFilePath(0, 0));
	if (referenceGenomeFile->is_open()) {
		*referenceGenomeFile >> nextReferenceGenomeSeq;
		referenceGenomeFile->close();
		delete referenceGenomeFile;
	}
	else {
		cout << "ReferenceGenome File �б� �Ұ�" << endl;
		exit(1);
	}

	for (int i = 0; i < FILE_SIZE; i++) {
		cout << i << "��° ReferenceGenome bit ��ȯ ��..!" << endl;
		string currentReferenceGenomeSeq = nextReferenceGenomeSeq;
		if (i != FILE_SIZE - 1) {
			referenceGenomeFile = new ifstream(createFilePath(0, i + 1));
			if (referenceGenomeFile->is_open()) {
				*referenceGenomeFile >> nextReferenceGenomeSeq;
				referenceGenomeFile->close();
				delete referenceGenomeFile;
			}
			else {
				cout << "ReferenceGenome File �б� �Ұ�" << endl;
				exit(1);
			}

			currentReferenceGenomeSeq += nextReferenceGenomeSeq.substr(0, L - 1);
		}

		unsigned long long int data = 0x0;
		for (int j = 0; j < currentReferenceGenomeSeq.length(); j++) {
			data <<= 2;
			data |= toBits(currentReferenceGenomeSeq[j]);
			if (j >= L - 1) {
				bitGenome[j - L + 1] = data;
			}
		}

		cout << i << "��° ReferenceGenome bit ��ȯ �Ϸ�..!" << endl;

		int matchedCount = 0;
		for (int j = 0; j < currentReferenceGenomeSeq.length() - L; j++) {
			set<unsigned long long int>::iterator it = shortRead.find(bitGenome[j]);
			if (it != shortRead.end()) {
				matchedCount = L;
				shortRead.erase(it);
			}

			if (matchedCount > 0) {
				matched[i][j] = true;
				matchedCount--;
			}
			else {
				mismatched[i].push_back(j);
			}
		}

		cout << i << "��° ReferenceGenome Matching check �Ϸ�..! mismatched rate " << mismatched[i].size() * 100 / N << "%" << endl;
	}
}

void checkMatchedRate() {
	int matchedCount[4] = { 0 };
	for (int i = 0; i < FILE_SIZE; i++) {
		string myGenomeSeq;
		ifstream *myGenomeFile = new ifstream(createFilePath(1, i));
		if (myGenomeFile->is_open()) {
			*myGenomeFile >> myGenomeSeq;
			myGenomeFile->close();
			delete myGenomeFile;
		}
		else {
			cout << "MyGenome File �б� �Ұ�" << endl;
			exit(1);
		}

		string referenceGenomeSeq;
		ifstream *referenceGenomeFile = new ifstream(createFilePath(0, i));
		if (referenceGenomeFile->is_open()) {
			*referenceGenomeFile >> referenceGenomeSeq;
			referenceGenomeFile->close();
			delete referenceGenomeFile;
		}
		else {
			cout << "ReferenceGenome File �б� �Ұ�" << endl;
			exit(1);
		}
		cout << i << "��° Matching check ��..!" << endl;

		for (int j = 0; j < N; j++) {
			if (matched[i][j] && referenceGenomeSeq[j] == myGenomeSeq[j]) {
				matchedCount[0]++;
			}
			else if (matched[i][j] && referenceGenomeSeq[j] != myGenomeSeq[j]) {
				matchedCount[1]++;
			}
			else if (!matched[i][j] && referenceGenomeSeq[j] == myGenomeSeq[j]) {
				matchedCount[2]++;
			}
			else if (!matched[i][j] && referenceGenomeSeq[j] != myGenomeSeq[j]) {
				matchedCount[3]++;
			}

		}
		cout << i << "��° Matching check �Ϸ�..!" << endl;
	}
	cout << "Matched ���¿��� ������ ������ ��� " << matchedCount[0] << endl;
	cout << "Matched ���¿��� ������ �ٸ� ��� " << matchedCount[1] << endl;
	cout << "Mismatched ���¿��� ������ ������ ��� " << matchedCount[2] << endl;
	cout << "Mismatched ���¿��� ������ �ٸ� ��� " << matchedCount[3] << endl;
}

int main() {
	chrono::system_clock::time_point start = chrono::system_clock::now();
	readShortRead();
	recoveryMatchedGenome();
	checkMatchedRate();
	chrono::system_clock::time_point end = chrono::system_clock::now();
	chrono::duration<double> sec = end - start;

	cout << endl << "N: " << TOTAL_LEN << ", M: " << M << ", L: " << L << endl;
	cout << "Binary Transformation Comparison �ɸ� �ð�: " << sec.count() << "��" << endl;
	return 0;
}