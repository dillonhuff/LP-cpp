#include <iostream>
#include <vector>

#include <gmpxx.h>

using namespace std;

struct value {

};

struct point {
  std::vector<value*> coords;
};

struct linear_expr {
  std::vector<value*> coeffs;
  value* constant;

  void set_coeff(const int dim, value* const v) {
    coeffs[dim] = v;
  }
};

struct basic_set {
  std::vector<linear_expr*> constraints;
};

struct context {

  value* val_alloc(const int val) {
    return new value();
  }

  linear_expr* linear_expr_alloc(const int dims) {
    auto e = new linear_expr();
    e->coeffs.resize(dims);
    return e;
  }

  basic_set* basic_set_alloc(const int dims) {
    return new basic_set();
  }
};

int main() {
  context ctx;

  auto a = ctx.basic_set_alloc(2);
  auto sum = ctx.linear_expr_alloc(2);
  sum->set_coeff(0, ctx.val_alloc(1));
  sum->set_coeff(1, ctx.val_alloc(1));
  cout << "Done" << endl;
}

