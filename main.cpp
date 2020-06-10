#include <cassert>
#include <iostream>
#include <vector>

#include <gmpxx.h>

using namespace std;

struct value {
  mpz_class v;

  value() : v(0) {}
  value(const int val) : v(val) {}
  value(const mpz_class& val) : v(val) {}
};

value* neg(value* v) {
  return new value(-(v->v));
}

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

  int dimension() const {
    return coeffs.size();
  }


  void set_const(value* const v) {
    constant = v;
  }

  value* get_const() {
    return constant;
  }

  value* get_coeff(const int dim) {
    assert(dim < dimension());
    return coeffs.at(dim);
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

struct standard_form {
  linear_expr* objective;
  vector<linear_expr*> equalities;
  int num_base_vars;
  int num_slack_vars;
};

standard_form to_standard_form(linear_expr* obj, const vector<linear_expr*>& constraints) {
  standard_form form;

  form.num_slack_vars = constraints.size();
  form.num_base_vars = 2*obj->dimension();
  int standard_dim = form.num_slack_vars + form.num_base_vars;

  form.objective = new linear_expr(standard_dim);
  for (int d = 0; d < obj->dimension(); d++) {
    form.objective->set_coeff(2*d, obj->get_coeff(d));
    form.objective->set_coeff(2*d + 1, neg(obj->get_coeff(d)));
  }

  int slack_offset = form.num_base_vars;
  for (auto c : constraints) {
    auto cs = new linear_expr(standard_dim);
    for (int d = 0; d < c->dimension(); d++) {
      cs->set_coeff(2*d, c->get_coeff(d));
      cs->set_coeff(2*d + 1, neg(c->get_coeff(d)));
    }

    cs->set_coeff(slack_offset, new value(1));

    cs->set_const(c->get_const());

    form.equalities.push_back(cs);
    slack_offset++;
  }

  return form;
}

value* maximize(linear_expr* sum, const vector<linear_expr*>& constraints) {
  cout << "Maximizing : " << *sum << endl;
  cout << "Subject to: " << endl;
  for (auto c : constraints) {
    cout << "  " << *c << " >= 0" << endl;
  }

  standard_form sf = to_standard_form(sum, constraints);
  cout << "Standard form..." << endl;
  cout << "  " << *(sf.objective) << endl;
  cout << "st" << endl;
  for (auto c : sf.equalities) {
    cout << "  " << *c << " = 0" << endl;
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

