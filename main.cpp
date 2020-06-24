#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include "algorithm.h"

#include <gmpxx.h>

using namespace dbhc;
using namespace std;

struct value {
  mpq_class v;

  value() : v(0) {}
  value(const int val) : v(val) {}
  value(const mpq_class& val) : v(val) {}
};

value* operator/(const value& v, const value& t) {
  return new value(v.v * t.v);
}
bool operator==(const value& v, const int t) {
  return cmp(v.v, t) == 0;
}

bool operator<(const value& v, const value& t) {
  return cmp(v.v, t.v) < 0;
}

bool operator>(const value& v, const value& t) {
  return cmp(v.v, t.v) > 0;
}

bool operator>(const value& v, const int t) {
  return cmp(v.v, t) > 0;
}

bool operator<(const value& v, const int t) {
  return cmp(v.v, t) < 0;
}

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

struct tableau {
  std::set<int> basic_variables;
  vector<vector<value*> > rows;

  value* constant(const int row) const {
    return get_entry(row, num_cols() - 1);
  }

  value* objective_coeff(const int col) const {
    return get_entry(0, col);
  }

  std::set<int> non_basic_variables() const {
    std::set<int> non_basic;
    for (int i = 0; i < num_cols() - 1; i++) {
      if (!elem(i, basic_variables)) {
        non_basic.insert(i);
      }
    }

    return non_basic;
  }

  tableau(const int nrows, const int ncols) {
    rows.resize(nrows);
    for (auto& r : rows) {
      r.resize(ncols);
    }
  }

  int num_rows() const {
    return rows.size();
  }

  int num_cols() const {
    assert(num_rows() > 0);
    return rows.at(0).size();
  }

  value* get_entry(const int r, const int c) const {
    assert(r < rows.size());
    assert(c < rows[r].size());
    return rows[r][c];
  }

  void set_entry(const int r, const int c, value* v) {
    assert(r < rows.size());
    assert(c < rows[r].size());
    rows[r][c] = v;
  }

};

tableau build_initial_tableau(standard_form& form) {
  int ncols = form.num_slack_vars + form.num_base_vars + 1;
  int nrows = form.equalities.size() + 1;
  tableau tab(nrows, ncols);

  for (int s = form.num_base_vars; s < ncols - 1; s++) {
    cout << s << " is a basic variable" << endl;
    tab.basic_variables.insert(s);
  }

  assert(tab.basic_variables.size() == form.equalities.size());

  int r = 0;
  for (int d = 0; d < form.objective->dimension(); d++) {
    tab.set_entry(r, d, form.objective->get_coeff(d));
  }
  tab.set_entry(r, ncols - 1, form.objective->get_const());

  r++;
  for (auto c : form.equalities) {
    for (int d = 0; d < c->dimension(); d++) {
      tab.set_entry(r, d, c->get_coeff(d));
    }
    tab.set_entry(r, ncols - 1, c->get_const());
    r++;
  }
  return tab;
}

int pick_pivot_row(tableau& tab) {
  int pr = 0;
  value* min = new value(0);
  for (int c = 0; c < tab.num_cols(); c++) {
    if (*(tab.get_entry(tab.num_rows() - 1, c)) > 0) {
      pr = c;
      min = tab.get_entry(tab.num_rows() - 1, c);
    }
  }
  return pr;
}

int pick_pivot_col(const int pivot_row, tableau& tab) {
  return 0;
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

  tableau tab = build_initial_tableau(sf);
  cout << "Tableau" << endl;
  cout << "  # basic vars: " << tab.basic_variables.size() << endl;
  for (auto b : tab.basic_variables) {
    cout << "    " << b << endl;
  }

  cout << "  # non basic vars: " << tab.non_basic_variables().size() << endl;
  for (auto b : tab.non_basic_variables()) {
    cout << "    " << b << endl;
    assert(!elem(b, tab.basic_variables));
  }
  for (int r = 0; r < tab.num_rows(); r++) {
    for (int c = 0; c < tab.num_cols(); c++) {
      cout << *(tab.get_entry(r, c)) << " ";
    }
    cout << endl;
  }

  int next_pivot_col = -1;
  value* old_ratio = new value(0);
  for (auto x : tab.non_basic_variables()) {
    value* c = tab.objective_coeff(x);
    cout << "a_" << x << " = " << *c << endl;
    if (*c > 0) {
      cout << "  " << "has positive coefficient" << endl;
      //int c = x;
      for (int r = 1; r < tab.num_rows(); r++) {
        cout << "   " << "b_" << r << " = " << *tab.constant(r) << endl;
        value* ratio_val = *(tab.constant(r)) / *c;
        cout << "    ratio = " << *ratio_val << endl;
        cout << "old ratio = " << *old_ratio << endl;
        if (*ratio_val > *old_ratio) {
          old_ratio = ratio_val;
          next_pivot_col = x;
        } else {
          cout << *ratio_val << " <= " << *old_ratio << endl;
        }
      }
    }
  }

  if (next_pivot_col == -1) {
    cout << "No pivot row found" << endl;
    assert(false);
  }


  cout << "next pivot col = " << next_pivot_col << endl;

  for (int r = 1; r < tab.num_rows(); r++) {
    value* b = tab.const_coeff(r);
    value* c = tab.variable_coeff(r, next_pivot_col);
    cout << "b = " << *b << endl;
    cout << "c = " << *c << endl;
  }

  //int pivot_row = pick_pivot_row(tab);
  //cout << "pivot row = " << pivot_row << endl;

  //int pivot_col = pick_pivot_col(pivot_row, tab);


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

