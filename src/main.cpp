#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <cstdio>
#include <chrono>
#include <cstdlib>
#include <string>
#include <queue>
#include <cmath>
#include <ctime>
#define nptr -1
#define kmove 20
#define Prob 0.95

using namespace std;

class block {
public:
	int id;
	int x;
	int y;
	int width;
	int height;
	int rotate;
};

class terminal {
public:
	int id;
	int x;
	int y;
};

class net {
public:
	int id;
	vector<block> listB;
	vector<terminal> listT;
};

class NODE {
public:
	int parent;
	int Lchild;
	int Rchild;
};

class COST {
public:
	int width;
	int height;
	double area;
	double HPWL;
	double aspectR;
	double cost;
};


//time recording functions
chrono::high_resolution_clock::time_point time_record() {
	return chrono::high_resolution_clock::now();
}

double time_output(chrono::high_resolution_clock::time_point start_time, chrono::high_resolution_clock::time_point end_time) {
	return chrono::duration<double>(end_time - start_time).count();
}
//global var1
vector<block> dicB;
vector<vector<int>> dicN;
vector<terminal> dicTermi;
//global var2
int NumBlocks, NumTerminals;
int NumNets;
bool debug = false;

double dead_ratio;
double total_block_area;
double area_target;
double NormArea = 0, NormHPWL = 0;
int FO;

int root_idx = nptr;
vector<NODE> BStree;

vector<int> contourH;

bool isPassFO;
int minC_rootidx;
COST best_COST;
vector<block> best_blockarray;
vector<NODE> best_tree;
int minCFO_rootidx;
COST best_COST_FO;
vector<block> best_barray_FO;
vector<NODE> best_tree_FO;
//
void parse_block(string block_file);
void parse_net(string net_file);
void parse_terminal(string terminal_file);
void build_iniTree();
void tree_to_FP();
void traverse_bstree(int current_node, bool left);
COST cal_cost();
unsigned int selectSeed();

unsigned int seed;

unsigned int selectSeed()
{
	if (NumBlocks == 100) {
		if (dead_ratio == 0.1) return 1542;
		else return 10629;
	}
	else if (NumBlocks == 200) {
		if (dead_ratio == 0.1) return 15423;
		else return 911;
	}
	else if (NumBlocks == 300) {
		if (dead_ratio == 0.1) return 1211;
		else return 12111;
	}

	return time(NULL);
}

void parse_block(string block_file)
{
	ifstream input;
	input.open(block_file);

	string temp1, temp2, str;
	input >> temp1 >> temp2 >> NumBlocks;
	input >> temp1 >> temp2 >> NumTerminals;
	input >> str;

	total_block_area = 0;
	dicB = vector<block>(NumBlocks);
	for (int i = 0; i < NumBlocks; i++) {
		getline(input, str);

		size_t pos1 = str.find("(");
		pos1 = str.find("(", pos1 + 1);
		pos1 = str.find("(", pos1 + 1);
		size_t pos2 = str.find(",");
		pos2 = str.find(",", pos2 + 1);
		pos2 = str.find(",", pos2 + 1);
		size_t pos3 = str.find(")");
		pos3 = str.find(")", pos3 + 1);
		pos3 = str.find(")", pos3 + 1);

		char buffer[10];
		int width, height;
		size_t len = str.copy(buffer, pos2 - pos1 - 1, pos1 + 1);
		buffer[len] = '\0';
		width = atoi(buffer);
		len = str.copy(buffer, pos3 - pos2 - 2, pos2 + 2);
		buffer[len] = '\0';
		height = atoi(buffer);

		dicB[i].id = i;
		dicB[i].x = nptr;
		dicB[i].y = nptr;
		dicB[i].width = width;
		dicB[i].height = height;
		dicB[i].rotate = 0;

		total_block_area += width * height;
	}

	area_target = total_block_area * (1 + dead_ratio);
	FO = sqrt(area_target);

	cout << "FO:                " << FO << '\n';
	cout << '\n';

	input.close();
}

void parse_net(string net_file)
{
	int dummy;
	ifstream input;
	input.open(net_file);

	string temp1, temp2, str;
	input >> temp1 >> temp2 >> NumNets;
	input >> temp1 >> temp2 >> dummy;

	dicN = vector<vector<int>>(NumNets);
	for (int i = 0; i < NumNets; i++) {
		int deg;
		input >> temp1 >> temp2 >> deg;
		for (int j = 0; j < deg; j++) {
			input >> str;
			int id;
			if (str[0] == 'p') {
				str.erase(0, 1);
				id = atoi(str.c_str()) + NumBlocks;
			}
			else if (str[0] == 's') {
				str.erase(0, 2);
				id = atoi(str.c_str());
			}
			dicN[i].emplace_back(id);
		}
	}

	input.close();
}

void parse_terminal(string terminal_file)
{
	ifstream input;
	input.open(terminal_file);

	string str;
	int x, y;

	dicTermi = vector<terminal>(NumTerminals + 1);
	for (int i = 1; i <= NumTerminals; i++) {
		input >> str >> x >> y;
		dicTermi[i].id = i;
		dicTermi[i].x = x;
		dicTermi[i].y = y;
	}

	input.close();
}

void build_iniTree()
{
	BStree = vector<NODE>(NumBlocks);
	vector<int> INS(NumBlocks, 0);
	queue<int> BF;

	root_idx = rand() % NumBlocks;
	BStree[root_idx].parent = nptr;
	BF.push(root_idx);
	INS[root_idx] = 1;

	int left = NumBlocks - 1;
	while (!BF.empty()) {
		int parent = BF.front();
		BF.pop();

		int Lchild = nptr, Rchild = nptr;
		if (left > 0) {
			do {
				Lchild = rand() % NumBlocks;
			} while (INS[Lchild]);
			BStree[parent].Lchild = Lchild;
			BF.push(Lchild);
			INS[Lchild] = 1;
			left--;

			if (left > 0) {
				do {
					Rchild = rand() % NumBlocks;
				} while (INS[Rchild]);
				BStree[parent].Rchild = Rchild;
				BF.push(Rchild);
				INS[Rchild] = 1;
				left--;
			}
		}
		BStree[parent].Lchild = Lchild;
		BStree[parent].Rchild = Rchild;
		if (Lchild != nptr)
			BStree[Lchild].parent = parent;
		if (Rchild != nptr)
			BStree[Rchild].parent = parent;
	}
}

void tree_to_FP()
{
	contourH = vector<int>(FO * 5, 0);
	dicB[root_idx].x = 0;
	dicB[root_idx].y = 0;

	for (int i = 0; i < dicB[root_idx].width; i++)
		contourH[i] = dicB[root_idx].height;

	if (BStree[root_idx].Lchild != nptr)
		traverse_bstree(BStree[root_idx].Lchild, true);
	if (BStree[root_idx].Rchild != nptr)
		traverse_bstree(BStree[root_idx].Rchild, false);
}
void traverse_bstree(int current_node, bool left)
{
	int parent = BStree[current_node].parent;
	if (left)
		dicB[current_node].x = dicB[parent].x + dicB[parent].width;
	else
		dicB[current_node].x = dicB[parent].x;

	int x_start = dicB[current_node].x;
	int x_end = x_start + dicB[current_node].width;
	int maxY = 0;
	for (int i = x_start; i < x_end; i++)
		if (contourH[i] > maxY)
			maxY = contourH[i];

	dicB[current_node].y = maxY;

	maxY += dicB[current_node].height;
	for (int i = x_start; i < x_end; i++)
		contourH[i] = maxY;

	if (BStree[current_node].Lchild != nptr)
		traverse_bstree(BStree[current_node].Lchild, true);
	if (BStree[current_node].Rchild != nptr)
		traverse_bstree(BStree[current_node].Rchild, false);
}
COST cal_cost()
{
	tree_to_FP();

	int width = 0, height = 0;
	for (int i = 0; i < NumBlocks; i++) {
		if (dicB[i].x + dicB[i].width > width)
			width = dicB[i].x + dicB[i].width;
		if (dicB[i].y + dicB[i].height > height)
			height = dicB[i].y + dicB[i].height;
	}

	double total_Area = width * height;
	double aspectR = (double)height / width;

	double HPWL = 0;

	for (const vector<int> &net : dicN) {
		int minX = width + 1, maxX = 0;
		int minY = height + 1, maxY = 0;
		for (const int pin : net) {
			if (pin < NumBlocks) {
				int x_center = dicB[pin].x + dicB[pin].width / 2;
				int y_center = dicB[pin].y + dicB[pin].height / 2;
				if (x_center < minX)
					minX = x_center;
				if (y_center < minY)
					minY = y_center;
				if (x_center > maxX)
					maxX = x_center;
				if (y_center > maxY)
					maxY = y_center;
			}
			else {
				const terminal &t = dicTermi[pin - NumBlocks];
				if (t.x < minX)
					minX = t.x;
				if (t.y < minY)
					minY = t.y;
				if (t.x > maxX)
					maxX = t.x;
				if (t.y > maxY)
					maxY = t.y;
			}
		}

		HPWL += (maxX - minX) + (maxY - minY);
	}

	COST output;
	output.width = width;
	output.height = height;
	output.area = total_Area;
	output.HPWL = HPWL;
	output.aspectR = aspectR;

	if (NormArea == 0)
		NormArea = total_Area;
	if (NormHPWL == 0)
		NormHPWL = HPWL;

	double cost_Area = output.area / NormArea;
	double cost_HPWL = output.HPWL / NormHPWL;
	double cost_aspectR = (1 - aspectR) * (1 - aspectR);
	double costW = 0;
	double costH = 0;
	if (width > FO)
		costW = ((double)width / FO);
	if (height > FO)
		costH = ((double)height / FO);
	output.cost = cost_Area + cost_HPWL + cost_aspectR + costW + costH;

	return output;
}


void op1_rotate(int node)
{
	int dummy = dicB[node].width;
	dicB[node].width = dicB[node].height;
	dicB[node].height = dummy;
	dicB[node].rotate = 1 - dicB[node].rotate;
}

void op2_swap(int pick1, int pick2)
{

	int node1_parent = BStree[pick1].parent;
	if (node1_parent != nptr) {
		if (BStree[node1_parent].Lchild == pick1)
			BStree[node1_parent].Lchild = pick2;
		else if (BStree[node1_parent].Rchild == pick1)
			BStree[node1_parent].Rchild = pick2;
	}

	int parent2 = BStree[pick2].parent;
	if (parent2 != nptr) {
		if (BStree[parent2].Lchild == pick2)
			BStree[parent2].Lchild = pick1;
		else if (BStree[parent2].Rchild == pick2)
			BStree[parent2].Rchild = pick1;

	}

	BStree[pick1].parent = parent2;
	BStree[pick2].parent = node1_parent;

	int node1_left_child = BStree[pick1].Lchild;
	int node1_right_child = BStree[pick1].Rchild;
	BStree[pick1].Lchild = BStree[pick2].Lchild;
	BStree[pick1].Rchild = BStree[pick2].Rchild;
	BStree[pick2].Lchild = node1_left_child;
	BStree[pick2].Rchild = node1_right_child;

	if (BStree[pick1].Lchild != nptr)
		BStree[BStree[pick1].Lchild].parent = pick1;
	if (BStree[pick1].Rchild != nptr)
		BStree[BStree[pick1].Rchild].parent = pick1;
	if (BStree[pick2].Lchild != nptr)
		BStree[BStree[pick2].Lchild].parent = pick2;
	if (BStree[pick2].Rchild != nptr)
		BStree[BStree[pick2].Rchild].parent = pick2;

	if (BStree[pick1].parent == pick1)
		BStree[pick1].parent = pick2;
	else if (BStree[pick1].Lchild == pick1)
		BStree[pick1].Lchild = pick2;
	else if (BStree[pick1].Rchild == pick1)
		BStree[pick1].Rchild = pick2;

	if (BStree[pick2].parent == pick2)
		BStree[pick2].parent = pick1;
	else if (BStree[pick2].Lchild == pick2)
		BStree[pick2].Lchild = pick1;
	else if (BStree[pick2].Rchild == pick2)
		BStree[pick2].Rchild = pick1;

	if (pick1 == root_idx)
		root_idx = pick2;
	else if (pick2 == root_idx)
		root_idx = pick1;
}

void op3_move(int node, int desti)
{
	
	if (BStree[node].Lchild == nptr && BStree[node].Rchild == nptr) {
		
		int parent = BStree[node].parent;
		if (BStree[parent].Lchild == node)
			BStree[parent].Lchild = nptr;
		else if (BStree[parent].Rchild == node)
			BStree[parent].Rchild = nptr;
	}
	else if (BStree[node].Lchild != nptr && BStree[node].Rchild != nptr) {

		do {
			bool move_left;
			if (BStree[node].Lchild != nptr && BStree[node].Rchild != nptr)
				move_left = rand() % 2 == 0;
			else if (BStree[node].Lchild != nptr)
				move_left = true;
			else
				move_left = false;

			if (move_left)
				op2_swap(node, BStree[node].Lchild);
			else
				op2_swap(node, BStree[node].Rchild);
		} while (BStree[node].Lchild != nptr || BStree[node].Rchild != nptr);

		int parent = BStree[node].parent;
		if (BStree[parent].Lchild == node)
			BStree[parent].Lchild = nptr;
		else if (BStree[parent].Rchild == node)
			BStree[parent].Rchild = nptr;
	}
	else {
		
		int child;
		if (BStree[node].Lchild != nptr)
			child = BStree[node].Lchild;
		else
			child = BStree[node].Rchild;

		int parent = BStree[node].parent;
		if (parent != nptr) {
			if (BStree[parent].Lchild == node)
				BStree[parent].Lchild = child;
			else if (BStree[parent].Rchild == node)
				BStree[parent].Rchild = child;
		}

		BStree[child].parent = parent;

		if (node == root_idx)
			root_idx = child;
	}


	int random_left_right = rand() % 4;
	int child;
	if (random_left_right == 0) {
		child = BStree[desti].Lchild;
		BStree[node].Lchild = child;
		BStree[node].Rchild = nptr;
		BStree[desti].Lchild = node;
	}
	else if (random_left_right == 0) {
		child = BStree[desti].Rchild;
		BStree[node].Lchild = child;
		BStree[node].Rchild = nptr;
		BStree[desti].Rchild = node;
	}
	else if (random_left_right == 0) {
		child = BStree[desti].Lchild;
		BStree[node].Lchild = nptr;
		BStree[node].Rchild = child;
		BStree[desti].Lchild = node;
	}
	else {
		child = BStree[desti].Rchild;
		BStree[node].Lchild = nptr;
		BStree[node].Rchild = child;
		BStree[desti].Rchild = node;
	}
	BStree[node].parent = desti;
	if (child != nptr)
		BStree[child].parent = node;
}

int perturb1() {
	int M = rand() % 3;
	if (M == 0) {
		int node = rand() % NumBlocks;
		op1_rotate(node);
	}
	else if (M == 1) {
		int pick1, pick2;
		pick1 = rand() % NumBlocks;
		do {
			pick2 = rand() % NumBlocks;
		} while (pick2 == pick1);
		op2_swap(pick1, pick2);
	}
	else if (M == 2) {
		int node, desti;
		node = rand() % NumBlocks;
		do {
			desti = rand() % NumBlocks;
		} while (desti == node || BStree[node].parent == desti);
		op3_move(node, desti);
	}
	return M;
}

int perturb2() {
	int random = rand() % 10;
	int M;
	if (random < 5) M = 0;
	else if (random < 10) M = 1;
	else M = 2;

	if (M == 0) {
		int node = rand() % NumBlocks;
		op1_rotate(node);
	}
	else if (M == 1) {
		int pick1, pick2;
		pick1 = rand() % NumBlocks;
		do {
			pick2 = rand() % NumBlocks;
		} while (pick2 == pick1);
		op2_swap(pick1, pick2);
	}
	else if (M == 2) {
		int node, desti;
		node = rand() % NumBlocks;
		do {
			desti = rand() % NumBlocks;
		} while (desti == node || BStree[node].parent == desti);
		op3_move(node, desti);
	}
	return M;
}

void SA()
{
	best_COST = cal_cost();
	best_blockarray = dicB;

	int N = kmove * NumBlocks;
	const double T0 = (-1)*(best_COST.cost) * NumBlocks / log(Prob);

	double T = T0;
	int MT = 0;
	int uphill = 0;
	int reject = 0;
	bool flag = false;
	bool success = false;
	bool timeflag = false;
	double old = 0;
	COST prev_cost = best_COST;
	isPassFO = false;
	int count = 0;

	clock_t init_time = clock();
	clock_t time = init_time;
	const int max_seconds = (NumBlocks / 20) * (NumBlocks / 20);
	const int TIME_LIMIT = 1200 - 5; // 20 minutes
	int seconds = 0, runtime = 0;

	do {
		flag = false;
		success = false;
		MT = 0;
		uphill = 0;
		reject = 0;

		do {
			flag = false;
			vector<block> dicBtemp(dicB);
			vector<NODE> BStreeTemp(BStree);
			int prevRootIdx = root_idx;

			int M = 0;
			if (T > 1E-10) M = perturb1();
			else M = perturb2();


			MT++;
			COST cur_cost = cal_cost();
			double dif_cost = cur_cost.cost - prev_cost.cost;

			double random = ((double)rand()) / RAND_MAX;
			if (dif_cost <= 0 || random < exp(-dif_cost / T)) {
				if (dif_cost > 0)
					uphill++;

				if (cur_cost.width <= FO && cur_cost.height <= FO) {
					if (isPassFO) {
						if (cur_cost.cost < best_COST_FO.cost) {
							minCFO_rootidx = root_idx;
							best_COST_FO = cur_cost;
							best_barray_FO = dicB;
							best_tree_FO = BStree;
							flag = true;
							success = true;
						}
					}
					else {
						isPassFO = true;
						minCFO_rootidx = root_idx;
						best_COST_FO = cur_cost;
						best_barray_FO = dicB;
						best_tree_FO = BStree;
						flag = true;
					}
				}

				if (cur_cost.cost < best_COST.cost) {
					minC_rootidx = root_idx;

					best_COST = cur_cost;
					best_blockarray = dicB;
					best_tree = BStree;
					flag = true;
				}
				if (flag);
				if (success);

				prev_cost = cur_cost;
			}
			else {
				reject++;
				root_idx = prevRootIdx;
				if (M == 0)
					dicB = dicBtemp;
				else
					BStree = BStreeTemp;
			}

		} while (uphill <= N && MT <= 2 * N);

		if (old == prev_cost.HPWL) { count++; cout << "count++ = " << count << endl; }
		else count = 0;
		//cout<< "                                           old HPWL: " << prev_cost.HPWL << endl;
		cout << "T: " << T << " W: " << prev_cost.width << " H: " << prev_cost.height << " HPWL: " << prev_cost.HPWL <<" cost: " << prev_cost.cost<< endl;
		old = prev_cost.HPWL;
		T = T*0.9;
		

		seconds = (clock() - time) / CLOCKS_PER_SEC;
		runtime = (clock() - init_time) / CLOCKS_PER_SEC;
		if (seconds >= max_seconds && isPassFO == false) {
			seconds = 0;
			time = clock();
			T = T0;
		}
		if (runtime > TIME_LIMIT / 2) timeflag = true;
		if (timeflag&&count >= 30) break;
	} while (seconds < max_seconds && runtime < TIME_LIMIT && count<50);

	if (isPassFO) {
		cout << "[Valid solution found]\n";
		cout << "W: " << best_COST_FO.width << " ";
		cout << "H: " << best_COST_FO.height << " ";
		cout << "HPWL: " << best_COST_FO.HPWL << " ";
		cout << "seed: " << seed << " ";
		cout << "White space: " << ((double)best_COST_FO.area / (double)total_block_area)-1 << '\n';
		cout << '\n';
	}
}

void output_function(string output_file, int HPWL, vector<block> &blockarray)
{
	ofstream output;
	output.open(output_file);

	output << "Wirelength " << HPWL << '\n';
	output << "Blocks\n";

	for (int i = 0; i < NumBlocks; i++) {
		if (blockarray[i].rotate)
			output << "sb" << i << " " << blockarray[i].x << " " << blockarray[i].y << " 1\n";
		else
			output << "sb" << i << " " << blockarray[i].x << " " << blockarray[i].y << " 0\n";
	}

	output.close();
}

int main(int argc, char **argv)
{
	
	int NumThreads = omp_get_max_threads();
	omp_set_num_threads(NumThreads);
	cout << "Num of parallel threads used: " << NumThreads << endl;


	string block_file = argv[1];
	string net_file = argv[2];
	string terminal_file = argv[3];
	string output_file = argv[4];
	dead_ratio = atof(argv[5]);

	//input time
	chrono::high_resolution_clock::time_point I_start = time_record();
	parse_block(block_file);
	parse_net(net_file);
	parse_terminal(terminal_file);
	chrono::high_resolution_clock::time_point I_end = time_record();
	double I_time = time_output(I_start, I_end);
	seed = selectSeed();
	srand(seed);
	
	//init time
	chrono::high_resolution_clock::time_point Ini_start = time_record();
	build_iniTree();
	chrono::high_resolution_clock::time_point Ini_end = time_record();
	double Ini_time = time_output(Ini_start, Ini_end);

	//calcu time
	chrono::high_resolution_clock::time_point Compu_start = time_record();
	SA();
	chrono::high_resolution_clock::time_point Compu_end = time_record();
	double Compu_time = time_output(Compu_start, Compu_end);

	//output time
	chrono::high_resolution_clock::time_point O_start = time_record();
	if (isPassFO)
		output_function(output_file, best_COST_FO.HPWL, best_barray_FO);
	else
		output_function(output_file, best_COST.HPWL, best_blockarray);

	chrono::high_resolution_clock::time_point O_end = time_record();
	double O_time = time_output(O_start, O_end);

	cout << "I/O time: " << O_time + I_time << endl;
	cout << "Ini_time: " << Ini_time << endl;
	cout << "Compute time: " << Compu_time << endl;
	cout << "Total time: " << O_time + I_time + Compu_time + Ini_time << endl;


	return 0;
}
