#include <iostream>
#include <vector>

#include <gmpxx.h>

using namespace std;

struct value {
  mpz_class v;

  value() : v(0) {}
  value(const int val) : v(val) {}
};

struct point {
  std::vector<value*> coords;
};

struct linear_expr {
  std::vector<value*> coeffs;
  value* constant;

  linear_expr(const int dimension) {
    coeffs.resize(dimension);
    for (int i = 0; i < dimension; i++) {
      coeffs[i] = new value(0);
    }
    constant = new value(0);
  }

  void set_const(value* const v) {
    constant = v;
  }

  void set_coeff(const int dim, value* const v) {
    coeffs[dim] = v;
  }
};

std::ostream& operator<<(std::ostream& out, const value& e) {
  out << e.v.get_str();
  return out;
}

std::ostream& operator<<(std::ostream& out, const linear_expr& e) {
  for (int i = 0; i < e.coeffs.size(); i++) {
    if (i > 0) {
      out << " + ";
    }
    out << *(e.coeffs.at(i)) << " x_" << i;
  }

  out << " + " << *(e.constant);
  return out;
}

struct basic_set {
  std::vector<linear_expr*> constraints;
};

struct context {

  value* val_alloc(const int val) {
    return new value(val);
  }

  linear_expr* linear_expr_alloc(const int dims) {
    auto e = new linear_expr(dims);
    return e;
  }

  basic_set* basic_set_alloc(const int dims) {
    return new basic_set();
  }
};

value* maximize(linear_expr* sum, const vector<linear_expr*>& constraints) {
  cout << "Maximizing : " << *sum << endl;
  cout << "Subject to: " << endl;
  for (auto c : constraints) {
    cout << "  " << *c << " >= 0" << endl;
  }
  return nullptr;
}

int main() {
  context ctx;

  auto a = ctx.basic_set_alloc(2);
  auto sum = ctx.linear_expr_alloc(1);
  sum->set_coeff(0, ctx.val_alloc(1));

  auto lc = ctx.linear_expr_alloc(1);
  lc->set_coeff(0, ctx.val_alloc(1));
  lc->set_const(ctx.val_alloc(1));

  vector<linear_expr*> constraints{lc};
  value* result = maximize(sum, constraints);

  cout << "Done" << endl;
}

