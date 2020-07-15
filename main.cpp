#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include "algorithm.h"

#include <gmpxx.h>

using namespace dbhc;
using namespace std;

template<typename T>
class matrix {

  vector<vector<T> > rows;

  void scale_row(const T& factor, const int row) {
    assert(row < num_rows());
    for (int c = 0; c < num_cols(); c++) {
      rows[row][c] = (rows[row][c]) * factor;
    }
  }

  T operator()(const int row, const int col) {
    return get_entry(row, col);
  }

  matrix(const int nrows, const int ncols) {
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

  T get_entry(const int r, const int c) const {
    assert(r < rows.size());
    assert(c < rows[r].size());
    return rows[r][c];
  }

  void set_entry(const int r, const int c, const T& v) {
    assert(r < rows.size());
    assert(c < rows[r].size());
    rows[r][c] = v;
  }
};

struct value {
  mpq_class v;

  value() : v(0) {}
  value(const int val) : v(val) {}
  value(const mpq_class& val) : v(val) {}
};

value operator-(const value& v, const value& t) {
  return value(v.v - t.v);
}

value operator+(const value& v, const value& t) {
  return value(v.v + t.v);
}

value operator*(const value& v, const value& t) {
  return value(v.v * t.v);
}

value operator/(const value& v, const value& t) {
  return value(v.v / t.v);
}

bool operator==(const value& v, const int t) {
  return cmp(v.v, t) == 0;
}

bool operator<(const value& v, const value& t) {
  return cmp(v.v, t.v) < 0;
}

bool operator<(const value& v, const int t) {
  return cmp(v.v, t) < 0;
}

bool operator>(const value& v, const value& t) {
  return cmp(v.v, t.v) > 0;
}

bool operator>(const value& v, const int t) {
  return cmp(v.v, t) > 0;
}

value neg(const value v) {
  return value(-(v.v));
}

struct point {
  std::vector<value*> coords;
};

struct linear_expr {
  std::vector<value> coeffs;
  value constant;

  linear_expr(const int dimension) {
    coeffs.resize(dimension);
    for (int i = 0; i < dimension; i++) {
      coeffs[i] = value(0);
    }
    constant = value(0);
  }

  int dimension() const {
    return coeffs.size();
  }


  void set_const(const value& v) {
    constant = v;
  }

  value get_const() {
    return constant;
  }

  value get_coeff(const int dim) {
    assert(dim < dimension());
    return coeffs.at(dim);
  }

  void set_coeff(const int dim, const value& v) {
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
    out << (e.coeffs.at(i)) << " x_" << i;
  }

  out << " + " << (e.constant);
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
  linear_expr objective;
  vector<linear_expr> equalities;
  int num_base_vars;
  int num_slack_vars;
};

standard_form to_standard_form(const linear_expr& obj, const vector<linear_expr>& constraints) {
  standard_form form;

  form.num_slack_vars = constraints.size();
  form.num_base_vars = 2*obj->dimension();
  int standard_dim = form.num_slack_vars + form.num_base_vars;

  form.objective = linear_expr(standard_dim);
  for (int d = 0; d < obj->dimension(); d++) {
    form.objective->set_coeff(2*d, obj->get_coeff(d));
    form.objective->set_coeff(2*d + 1, neg(obj->get_coeff(d)));
  }

  int slack_offset = form.num_base_vars;
  for (auto c : constraints) {
    auto cs = linear_expr(standard_dim);
    for (int d = 0; d < c->dimension(); d++) {
      cs->set_coeff(2*d, c->get_coeff(d));
      cs->set_coeff(2*d + 1, neg(c->get_coeff(d)));
    }

    cs->set_coeff(slack_offset, value(1));

    cs->set_const(c->get_const());

    form.equalities.push_back(cs);
    slack_offset++;
  }

  return form;
}

struct tableau {
  value maximum;
  std::set<int> basic_variables;
  vector<vector<value> > rows;

  void subtract_row(const std::vector<value>& diffs, const int row) {
    assert(diffs.size() == num_cols());
    for (int c = 0; c < num_cols(); c++) {
      set_entry(row, c, get_entry(row, c) - diffs.at(c));
    }
  }

  void scale_row(const value& factor, const int row) {
    assert(row < num_rows());
    for (int c = 0; c < num_cols(); c++) {
      rows[row][c] = (rows[row][c]) * factor;
    }
  }

  void exchange(const int new_basic, const int old_basic) {
    assert(elem(old_basic, basic_variables));
    basic_variables.insert(new_basic);
    basic_variables.erase(old_basic);
  }

  value const_coeff(const int row) const {
    return constant(row);
  }

  value constant(const int row) const {
    return get_entry(row, num_cols() - 1);
  }

  value variable_coeff(const int row, const int col) const {
    return get_entry(row, col);
  }

  value objective_coeff(const int col) const {
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
    maximum = value(0);
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

  value get_entry(const int r, const int c) const {
    assert(r < rows.size());
    assert(c < rows[r].size());
    return value(rows[r][c].v);
  }

  void set_entry(const int r, const int c, const value& v) {
    assert(r < rows.size());
    assert(c < rows[r].size());
    rows[r][c] = v;
  }

  void print(std::ostream& out) {
    for (int r = 0; r < num_rows(); r++) {
      for (int c = 0; c < num_cols(); c++) {
        cout << get_entry(r, c) << " ";
      }
      cout << endl;
    }
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
  for (int d = 0; d < form.objective.dimension(); d++) {
    tab.set_entry(r, d, form.objective.get_coeff(d));
  }
  tab.set_entry(r, ncols - 1, form.objective.get_const());

  r++;
  for (auto c : form.equalities) {
    for (int d = 0; d < c.dimension(); d++) {
      tab.set_entry(r, d, c.get_coeff(d));
    }
    tab.set_entry(r, ncols - 1, c.get_const());
    r++;
  }
  return tab;
}

int pick_pivot_row(const int next_pivot_col, tableau& tab) {
  value max(0);
  int pivot_row = -1;
  for (int r = 1; r < tab.num_rows(); r++) {
    value b = tab.const_coeff(r);
    value c = tab.variable_coeff(r, next_pivot_col);
    cout << "b = " << b << endl;
    cout << "c = " << c << endl;
    if (pivot_row == -1 || (b / c) > max) {
      max = (b / c);
      pivot_row = r;
    }
  }
  return pivot_row;
}

int pick_pivot_col(tableau& tab) {
  int next_pivot_col = -1;

  value min_obj_coeff(-1);

  for (int x = 0; x < tab.num_cols() - 1; x++) {
    value c = tab.objective_coeff(x);
    if (c > 0) {
      if (min_obj_coeff < 0 || c < min_obj_coeff) {
        min_obj_coeff = c;
        next_pivot_col = x;
      }
    }
  }

  if (next_pivot_col == -1) {
    cout << "No pivot row found" << endl;
    assert(false);
  }

  return next_pivot_col;
}

bool can_improve(tableau& tab) {
  for (int c = 0; c < tab.num_cols() - 1; c++) {
    if (tab.objective_coeff(c) > (int) 0) {
      return true;
    }
  }
  return false;
}

value maximize(linear_expr sum, const vector<linear_expr>& constraints) {
  cout << "Maximizing : " << *sum << endl;
  cout << "Subject to: " << endl;
  for (auto c : constraints) {
    cout << "  " << c << " >= 0" << endl;
  }

  standard_form sf = to_standard_form(sum, constraints);
  cout << "Standard form..." << endl;
  cout << "  " << (sf.objective) << endl;
  cout << "st" << endl;
  for (auto c : sf.equalities) {
    cout << "  " << c << " = 0" << endl;
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
      cout << (tab.get_entry(r, c)) << " ";
    }
    cout << endl;
  }

  while (can_improve(tab)) {
    int next_pivot_col = pick_pivot_col(tab);
    cout << "next pivot col = " << next_pivot_col << endl;

    int pivot_row = pick_pivot_row(next_pivot_col, tab);
    cout << "Pivot row = " << pivot_row << endl;

    assert(pivot_row >= 0);

    cout << "--- Next non basic variable to convert: " << next_pivot_col << endl;
    int next_pivot_row = -1;
    for (int r = 1; r < tab.num_rows(); r++) {
      next_pivot_row = r;
    }

    auto pivot_val = tab.get_entry(next_pivot_row, next_pivot_col);
    tab.maximum = tab.maximum - ((tab.objective_coeff(next_pivot_col)) * ((tab.const_coeff(next_pivot_row)) / pivot_val));
    cout << "Pivot val = " << pivot_val << endl;
    cout << "New max   = " << (tab.maximum) << endl;
    tab.scale_row((value(1) / pivot_val), pivot_row);

    for (int r = 0; r < tab.num_rows(); r++) {
      if (r != next_pivot_row) {
        vector<value> muls;
        for (int c = 0; c < tab.num_cols(); c++) {
          muls.push_back(tab.get_entry(r, c) * tab.get_entry(pivot_row, c));
        }
        tab.subtract_row(muls, r);
      }
    }

    tab.print(cout);
  }

  cout << "Final tableau" << endl;
  tab.print(cout);
  return tab.maximum;
}

void basic_test() {
  context ctx;

  linear_expr sum(1);
  sum.set_coeff(0, value(1));

  linear_expr lc(1);
  lc.set_coeff(0, value(-1));
  lc.set_const(value(5));

  vector<linear_expr> constraints{lc};
  value result = maximize(sum, constraints);

  assert(result == 5);

  cout << "BASIC TEST Passed" << endl;
}

void ft_test() {

}

int main() {
  basic_test();
  ft_test();
}

