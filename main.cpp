#include <iostream>
#include <vector>

using namespace std;

struct value {

};

struct linear_expr {
  std::vector<value*> coeffs;
  value* constant;
};

int main() {
  cout << "Done" << endl;
}
