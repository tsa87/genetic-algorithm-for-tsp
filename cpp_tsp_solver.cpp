#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <bitset>

using namespace std;
#define PROBLEM_SIZE 8

//compute the euclidian distance bewteen 2 points
float compute_dist(float x1, float y1, float x2, float y2) {
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

//return a n * n distance matrix
float** load_tsp(char* file_path, int n){

  int index;
  float x, y;
  ifstream inFile(file_path);

  //malloc
  float** cord_mat = new float*[n];
  float** dist_mat = new float*[n];
  for(int i = 0; i < n; i++){
    dist_mat[i] = new float[n];
  }

  while (inFile >> index >> x >> y) {
    cord_mat[index-1] = new float[2];
    cord_mat[index-1][0] = x;
    cord_mat[index-1][1] = y;
  }

  for(int i = 0; i < n; i++){
    for (int j = i; j < n; j++) {
      float distance = compute_dist(cord_mat[i][0], cord_mat[i][1], cord_mat[j][0], cord_mat[j][1]);
      dist_mat[i][j] = distance;
      dist_mat[j][i] = distance;
    }
  }

  for(int i = 0; i < n; i++) delete[] cord_mat[i];
  delete[] cord_mat;

  return dist_mat;
}

// dynamically generates the next subset representation
// ordered by the n-ary number created by the concatenation of sorted element of a subset
// ex. choosing 2 elements from {0, 1, 2, 3}
// set order: {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
// binary eqv: 1100, 1010, 1001, 0110, 0101, 0011
string increment(string& binary) {

  int size = binary.length();
  int i = size-1;
  int one_count;

  //Moves a 1 one spot forward.
  while (i >= 0) {
    if (binary[i] =='1') {
      one_count++;
      if (i != size-1){
        if (binary[i+1] == '0'){
          binary[i+1] = '1';
          binary[i] = '0';
          break;
        }
      }
    }
    i--;
  }

  //reorder permutation to be the next smallest binary string larger than input.
  i = i + 1; //move to the most recent 1
  int one_end_index = i + one_count -1;
  for (; i < size; i++) {
    if (i <= one_end_index) binary[i] = '1';
    else binary[i] = '0';
  }

  return binary;
}

void display_cost(float** cost_mat, int irange, int jrange){
  for (int i = 0; i < irange; i ++) {
    for (int j = 0; j < jrange; j++) {
      cout << cost_mat[i][j] << " ";
    }
    cout << endl << endl << endl;
  }
}


int main() {
  char file_path[] = "cities8.txt";
  int n = 9;

  // permutation is an interger whose binary representation indicates whether an element is part of a set
  // Ex. 9 -> 1001 -> {0, 3} is our set
  // NOTE: The starting point is always in our set implicitly.
  unsigned long long int permutation = (1 << (n-1)); //[0, 2^n-1]

  float** dist_mat = load_tsp(file_path, n);

  //initalize the 2^n * n-1 matrix with inf
  float** cost_mat = new float*[permutation];
  for (unsigned long long int i = 0; i < permutation; i++){
    cost_mat[i] = new float[n-1];
  }
  for (unsigned long long int i = 1; i < permutation; i++){
    for (int j = 0; j < n; j++) {
      cost_mat[i][j] = numeric_limits<float>::max();
    }
  }

  //backtrace matrix
  float** parent_mat = new float*[permutation];
  for (unsigned long long int i = 0; i < permutation; i++){
    parent_mat[i] = new float[n-1];
  }
  for (unsigned long long int i = 1; i < permutation; i++){
    for (int j = 0; j < n; j++) {
      parent_mat[i][j] = -1;
    }
  }

  // initalize distance of size 1 sets.
  // NOTE: Index 0 represents the second node (the one after starting node) for permutation purposes
  // NOTE: (with the exception of distance matrix)
  for (int i = 0; i < n-1; i++){
    cost_mat[1<<i][i] = dist_mat[0][i+1];
  }

  for (int size = 2; size < n; size++) { //Size not including Node 0

    //subset initalized to 11...100...0 where |1| = k and len(string) = n-1
    string subset = bitset<PROBLEM_SIZE>(((1 << size) - 1) << (n - 1 - size)).to_string();
    unsigned long long int permutation_id;

    while (permutation_id != ((1 << size) - 1)) {
      permutation_id = stoull(subset, nullptr, 2);

      //if j E S
      for (int j = 0; j < n-1; j++) {
        if (subset[j] == '1') {

          //generate the subset of S - {j} id
          string without_j = subset;
          without_j[j] = '0';
          unsigned long long int permutation_id_without_j = stoull(without_j, nullptr, 2);

          float min_cost = numeric_limits<float>::max();
          float cost;
          int best_i;

          // if i E S
          for (int i = 0; i < n-1; i++) {
            if (subset[i] == '1') {
              cost = cost_mat[permutation_id_without_j][i] + dist_mat[j][i];
              if (cost < min_cost) {
                min_cost = cost;
                best_i = i;
              }
            }
          }
          cost_mat[permutation_id][j] = min_cost;
          parent_mat[permutation_id][j] = best_i;
        }
      }
      subset = increment(subset);
    }
  }

  float min_cost = numeric_limits<float>::max();
  unsigned long long int permutation_id = permutation-1;
  float cost;
  int stop;
  for (int i = 0; i < n; i++) {
    cost = cost_mat[permutation_id][i] + dist_mat[0][i+1];
    if (cost < min_cost) {
      min_cost = cost;
      stop = i;
    }
  }

  //display_cost(cost_mat, permutation, n);
  cout << "Route Cost: " << min_cost<< endl;
  cout << "Route: " << "START" << " -> ";

  min_cost = numeric_limits<float>::max();
  string subset;
  while (true){
    stop = parent_mat[permutation_id][stop];
    if (stop == -1) break;

    cout << stop << " -> ";

    subset = bitset<PROBLEM_SIZE>(permutation_id).to_string();
    subset[stop] = '0';
    permutation_id = stoull(subset, nullptr, 2);

    for (int i = 0; i < n; i++) {
      cost = cost_mat[permutation_id][i];
      if (cost < min_cost) {
        min_cost = cost;
        stop = i;
      }
    }
  }

  cout << "START" << endl;

}
