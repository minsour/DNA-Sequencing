///* 2013113506 박민수 */
//#include <iostream>
//#include <sstream>
//#include <string>
//#include <vector>
//#include <random>
//#include <fstream>
//#include <list>
//#include <map>
//using namespace std;
//#define FILE_SIZE 50
//#define N 2000000
//#define M 25000000
//#define TOTAL_LEN N*FILE_SIZE
//#define DIFF_PERCENT 0.01
//#define L 32
//
/* 2013113506 박민수 */
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
//#define FILE_SIZE 50
//#define N 2000000
//#define M 25000000
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
set<unsigned long long int> shortRead;
unsigned long long int bitGenome[N];
vector<int> mismatched[FILE_SIZE];
bool matched[FILE_SIZE][N];
int shortReadIndex[M];

string createFilePath(int fileType, int fileIndex) {
	return folderName[fileType] + "\\genome" + to_string(fileIndex) + ".txt";
}

int getRandomNumber(int min, int max) {
	//< 1단계. 시드 설정
	random_device rn;
	mt19937_64 rnd(rn());

	//< 2단계. 분포 설정 ( 정수 )
	uniform_int_distribution<int> range(min, max);

	//< 3단계. 값 추출
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
		cout << "길이 " << N << "인 ReferenceGenomeSeq" << i << " 생성 중..!" << endl;
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
			cout << "ReferenceGenome File 생성 불가" << endl;
			exit(1);
		}
		cout << "ReferenceGenomeSeq" << i << " 생성 완료!" << endl;
	}
}

void makeMyGenome() {
	cout << "MyGenome 변경 Index 생성 중..!" << endl;
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
	cout << "MyGenome 변경 Index 생성 완료..!" << endl;

	int index = 0;
	for (int i = 0; i < FILE_SIZE && index < diffSize; i++) {
		cout << i << "번째 ReferenceGenome의 MyGenome 생성 중..!" << endl;

		ifstream *referenceGenomeFile = new ifstream(createFilePath(0, i));
		string genome;
		if (referenceGenomeFile->is_open()) {
			*referenceGenomeFile >> genome;
			referenceGenomeFile->close();
			delete referenceGenomeFile;
		}
		else {
			cout << "ReferenceGenome File 읽기 불가" << endl;
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
				cout << "MyGenome File 생성 불가" << endl;
				exit(1);
			}
			cout << i << "번째 ReferenceGenome의 MyGenome 생성 완료..!" << endl;
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
			cout << "MyGenome File 생성 불가" << endl;
			exit(1);
		}

		cout << i << "번째 ReferenceGenome의 MyGenome 생성 완료..!" << endl;
	}
}

void makeShortReadGenome() {
	cout << "ShortReadGenome Index 생성 중..!" << endl;

	for (int i = 0; i < M; i++) {
		shortReadIndex[i] = getRandomNumber(0, TOTAL_LEN - L);
	}

	sort(shortReadIndex, shortReadIndex + M);
	cout << "ShortReadGenome Index 생성 완료..!" << endl;

	int index = 0;
	ofstream *shortReadFile = new ofstream(createFilePath(2, 0));
	if (shortReadFile->is_open() == false) {
		cout << "Short Read File 생성 불가" << endl;
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
			cout << "ReferenceGenome File 읽기 불가" << endl;
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
					cout << "ReferenceGenome File 읽기 불가" << endl;
					exit(1);
				}
				genome += temp;
			}
			*shortReadFile << genome.substr(shortReadIndex[index] - i * N, L) << endl;
			index++;
		}
	}
	shortReadFile->close();
	delete shortReadFile;

	shortReadFile = new ofstream(createFilePath(2, index));
	for (int i = 0; i < M; i++) {
		*shortReadFile << shortReadIndex[i] << endl;
	}
	shortReadFile->close();
	delete shortReadFile;
}


void readShortRead() {
	ifstream *shortReadGenomeFile = new ifstream(createFilePath(2, 0));
	if (shortReadGenomeFile->is_open()) {
		for (int i = 0; i < M; i++) {
			if (i % 1000000 == 0) {
				cout << i << "/" << M << " ShortRead 읽는 중..." << endl;
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
		cout << "ShortRead Genome File 읽기 불가" << endl;
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
		cout << "ReferenceGenome File 읽기 불가" << endl;
		exit(1);
	}

	for (int i = 0; i < FILE_SIZE; i++) {
		cout << i << "번째 ReferenceGenome bit 변환 중..!" << endl;
		string currentReferenceGenomeSeq = nextReferenceGenomeSeq;

		if (i != FILE_SIZE - 1) {
			referenceGenomeFile = new ifstream(createFilePath(0, i + 1));
			if (referenceGenomeFile->is_open()) {
				*referenceGenomeFile >> nextReferenceGenomeSeq;
				referenceGenomeFile->close();
				delete referenceGenomeFile;
			}
			else {
				cout << "ReferenceGenome File 읽기 불가" << endl;
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

		cout << i << "번째 ReferenceGenome bit 변환 완료..!" << endl;

		int matchedCount = 0;
		for (int j = 0; j < currentReferenceGenomeSeq.length() - L; j++) {
			set<unsigned long long int>::iterator it = shortRead.find(bitGenome[j]);
			if (it != shortRead.end()) {
				matchedCount = L;
			}

			if (matchedCount > 0) {
				matched[i][j] = true;
				matchedCount--;
			}
			else {
				mismatched[i].push_back(j);
			}
		}

		cout << i << "번째 ReferenceGenome Matching check 완료..! mismatched rate " << mismatched[i].size() * 100 / N << "%" << endl;
	}
}


void printResult(vector<int> result, const string &text, const string &pattern) {
	for (int r : result) {
		cout << "Text String " << r << "번째에서 " << r + pattern.size() - 1
			<< "번째: " << text.substr(r, pattern.size()) << endl;
	}
}

int rabinKarp(string text, string pattern, vector<int> &result) {
	int d = 26, q = 13, p = 0, t = 0, patternLen = pattern.size(), textLen = text.size();
	int dh = (int)pow(d, patternLen - 1) % q;

	for (int i = 0; i <= patternLen - 1; i++) {
		p = (d * p + pattern[i]) % q;
		t = (d * t + text[i]) % q;
	}

	for (int s = 0; s < textLen - patternLen + 1; s++) {
		if (p == t) {
			int j;
			for (j = 0; j < patternLen; j++) {
				if (pattern[j] != text[s + j])
					break;
			}
			if (j == patternLen)
				result.push_back(s);
		}
		else if (s < textLen - patternLen)
			t = (d * (t - text[s] * dh) + text[s + patternLen]) % q;
	}

	return result.size();
}

void runRabinKarp() {
	string reference;
	ifstream *shortReadGenomeFile = new ifstream(createFilePath(2, 0));
	if (shortReadGenomeFile->is_open()) {
	}
	else {
		cout << "ReferenceGenome File 읽기 불가" << endl;
		exit(1);
	}
	ifstream *referenceGenomeFile = new ifstream(createFilePath(0, 0));
	if (referenceGenomeFile->is_open()) {
		*referenceGenomeFile >> reference;
		referenceGenomeFile->close();
		delete referenceGenomeFile;
	}
	else {
		cout << "ReferenceGenome File 읽기 불가" << endl;
		exit(1);
	}
	int sum = 0;
	for (int i = 0; i < M; i++) {
		string shortRead;
		*shortReadGenomeFile >> shortRead;
		string textStr, patternStr;
		vector<int> rabinResult;

		textStr = reference;
		patternStr = shortRead;
		/*cout << "Text String을 입력하세요: ";
		getline(cin, textStr);
		cout << "Pattern String을 입력하세요: ";
		getline(cin, patternStr);
		cout << endl;*/

		//chrono::system_clock::time_point start = chrono::system_clock::now();
		int result = rabinKarp(textStr, patternStr, rabinResult);
		//chrono::system_clock::time_point end = chrono::system_clock::now();
		//chrono::duration<double> sec = end - start;

		//printResult(rabinResult, textStr, patternStr);
		sum += result;
		if (i % 100 == 0) cout << i << ": 일치된 개수:" << sum << endl;
		//cout << "일치된 String 개수: " << result << endl;
		//cout << "Rabin Karp 걸린 시간: " << sec.count() << "초" << endl << endl;
	}

	shortReadGenomeFile->close();
	delete shortReadGenomeFile;
}

int main() {
	makeReferenceGenome();
	makeMyGenome();
	makeShortReadGenome();
	/*readShortRead();
	recoveryMatchedGenome();
	recoveryMismatchedGenome();*/

	chrono::system_clock::time_point start = chrono::system_clock::now();
	runRabinKarp();
	chrono::system_clock::time_point end = chrono::system_clock::now();
	chrono::duration<double> sec = end - start;

	cout << "N: " << N << ", M: " << M << ", L: " << L << endl;
	cout << "Rabin Karp 걸린 시간: " << sec.count() << "초" << endl << endl;

	return 0;
}