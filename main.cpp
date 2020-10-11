#include <cassert>
#include <iostream>
#include <iomanip>
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

bool operator==(const value& v, const int t) {
  return cmp(v.v, t) == 0;
}

bool operator!=(const value& v, const int t) {
  return !(v == t);
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

bool operator>=(const value& v, const int t) {
  return cmp(v.v, t) >= 0;
}

value operator-(const value& v, const value& t) {
  return value(v.v - t.v);
}

value operator+(const value& v, const value& t) {
  return value(v.v + t.v);
}

//value operator+=(const value& v, const value& t) {
  //return v + t;
//}

value operator*(const value& v, const value& t) {
  return value(v.v * t.v);
}

value operator/(const value& v, const value& t) {
  assert(t != 0);
  return value(v.v / t.v);
}

value neg(const value v) {
  return value(-(v.v));
}

struct point {
  std::vector<value> coords;

  point(const std::vector<value>& coords_) : coords(coords_) {}
  point(const value& a, const value& b) : coords({a, b}) {}
  point(const value& a) : coords({a}) {}

  const value& operator[](const int i) const {
    assert(i < dimension());
    return coords[i];
  }

  int dimension() const { return coords.size(); }
};

struct linear_expr {
  std::vector<value> coeffs;
  value constant_val;

  linear_expr() {
    coeffs.resize(0);
    constant_val = value(0);
  }

  linear_expr(const int dimension) {
    coeffs.resize(dimension);
    for (int i = 0; i < dimension; i++) {
      coeffs[i] = value(0);
    }
    constant_val = value(0);
  }

  int dimension() const {
    return coeffs.size();
  }

  const value& operator[](const int i) const {
    assert(i < dimension());
    return coeffs[i];
  }

  linear_expr scale(const int v) {
    linear_expr scaled = *this;
    for (int d = 0; d < dimension(); d++) {
      scaled.set_coeff(d, scaled.get_coeff(d)*v);
    }
    scaled.set_const(scaled.get_const()*v);
    return scaled;
  }

  void set_const(const value& v) {
    constant_val = v;
  }

  value constant() const {
    return constant_val;
  }

  value get_const() const {
    return constant();
  }

  value get_coeff(const int dim) const {
    assert(dim < dimension());
    return coeffs.at(dim);
  }

  void set_coeff(const int dim, const value& v) {
    coeffs[dim] = v;
  }
};

enum linear_constraint_type {
  LINEAR_CONSTRAINT_TYPE_LEQ,
  LINEAR_CONSTRAINT_TYPE_GEQ,
  LINEAR_CONSTRAINT_TYPE_EQ
};

struct linear_constraint {
  linear_constraint_type tp;
  linear_expr expr;

  value constant() const { return expr.get_const(); }

  int dimension() const { return expr.dimension(); }

  const value& operator[](const int i) const {
    assert(i < dimension());
    return expr[i];
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

  out << " + " << (e.constant());
  return out;
}

std::ostream& operator<<(std::ostream& out, const linear_constraint& lc) {
  out << lc.expr;
  if (lc.tp == LINEAR_CONSTRAINT_TYPE_EQ) {
    out << " = ";
  } else if (lc.tp == LINEAR_CONSTRAINT_TYPE_LEQ) {
    out << " <= ";
  } else if (lc.tp == LINEAR_CONSTRAINT_TYPE_GEQ) {
    cout << " >= ";
  }
  out << "0";
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

void sanity_check(standard_form& sf) {
  for (auto& e : sf.equalities) {
    assert(e.get_const() < 0);
  }
}

standard_form to_standard_form(const linear_expr& obj, const vector<linear_constraint>& constraints) {
  standard_form form;

  form.num_slack_vars = constraints.size();
  form.num_base_vars = 2*obj.dimension();
  int standard_dim = form.num_slack_vars + form.num_base_vars;

  form.objective = linear_expr(standard_dim);
  for (int d = 0; d < obj.dimension(); d++) {
    form.objective.set_coeff(2*d, obj.get_coeff(d));
    form.objective.set_coeff(2*d + 1, neg(obj.get_coeff(d)));
  }

  int slack_offset = form.num_base_vars;
  for (auto c : constraints) {
    auto cs = linear_expr(standard_dim);
    for (int d = 0; d < c.expr.dimension(); d++) {
      cs.set_coeff(2*d, c.expr.get_coeff(d));
      cs.set_coeff(2*d + 1, neg(c.expr.get_coeff(d)));
    }

    cs.set_coeff(slack_offset, value(1));

    cs.set_const(c.expr.get_const());
    //linear_constraint_type tp = c.tp;
    if (cs.get_const() >= 0) {
      ////cout << "Error: Standard form requires the constant to be negative when on LHS of the equality, but we have: " << cs << endl;
      ////assert(false);
      cs = cs.scale(-1);
    }

    form.equalities.push_back(cs);
    slack_offset++;
  }

  return form;
}

struct tableau {
  value maximum;
  std::set<int> basic_variables;
  vector<vector<value> > rows;

  value& operator()(const int r, const int c) {
    return rows[r][c];
  }

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
    //cout << "adding " << new_basic << " to basis, and removing: " << old_basic << endl;
    //assert(elem(old_basic, basic_variables));
    //basic_variables.insert(new_basic);
    //basic_variables.erase(old_basic);
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
    out << "Basic Variables: ";
    for (auto b : basic_variables) {
      out << b << " ";
    }
    out << endl;
    for (int r = 0; r < num_rows(); r++) {
      if (r == 0) {
        out << "--- Obj" << endl;
      }
      for (int c = 0; c < num_cols(); c++) {
        out << get_entry(r, c) << " ";
      }
      out << endl;
      if (r == 0) {
        out << "--- Constraints" << endl;
      }
    }
  }

};

std::ostream& operator<<(std::ostream& out, const tableau& tab) {
  int i = 0;
  for (auto& r : tab.rows) {
    for (auto& v : r) {
      std::cout << std::setfill(' ') << std::setw(8) << v << " ";
    }
    cout << endl;
    if (i == 0) {
      std::cout << std::setfill('-') << std::setw(8*(r.size() + 1)) << "" << std::endl;
    }
    i++;
  }
  return out;
}

void sanity_check(tableau& tab) {
  for (int r = 0; r < tab.num_rows(); r++) {
    assert(tab.const_coeff(r) >= 0);
  }
}

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
    tab.set_entry(r, ncols - 1, -1*c.get_const());
    r++;
  }
  return tab;
}

int pick_pivot_row(const int next_pivot_col, tableau& tab) {
  cout << "picking pivot row for " << next_pivot_col << endl;
  value max(0);
  int pivot_row = -1;
  for (int r = 1; r < tab.num_rows(); r++) {
    value b = tab.const_coeff(r);
    value c = tab.variable_coeff(r, next_pivot_col);
    if (c > 0) {
      cout << "b = " << b << endl;
      cout << "c = " << c << endl;
      if (pivot_row == -1 || (b / c) > max) {
        max = (b / c);
        pivot_row = r;
      }
    }
  }
  return pivot_row;
}

int pick_pivot_col(tableau& tab) {
  int next_pivot_col = -1;

  value min_obj_coeff(-1);

  for (int x = 0; x < tab.num_cols() - 1; x++) {
    value c = tab.objective_coeff(x);
    if (c < 0) {
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
    if (tab.objective_coeff(c) < (int) 0) {
      return true;
    }
  }
  return false;
}

enum lp_result_type {
  LP_RESULT_TYPE_OPTIMAL,
  LP_RESULT_TYPE_UNBOUNDED,
  LP_RESULT_TYPE_INFEASIBLE
};

struct lp_result {
  value val;
  lp_result_type tp;
};

lp_result maximize(linear_expr sum, const vector<linear_constraint>& constraints) {

  cout << "Maximizing : " << sum << endl;
  cout << "Subject to: " << endl;
  for (auto c : constraints) {
    cout << "  " << c << endl;
  }

  standard_form sf = to_standard_form(sum, constraints);
  cout << "Standard form..." << endl;
  cout << "  " << (sf.objective) << endl;
  cout << "st" << endl;
  for (auto c : sf.equalities) {
    cout << "  " << c << " = 0" << endl;
  }
  sanity_check(sf);

  tableau tab = build_initial_tableau(sf);

  cout << "Tableau" << endl;
  cout << tab << endl;
  assert(false);

  //cout << "  # non basic vars: " << tab.non_basic_variables().size() << endl;
  //for (auto b : tab.non_basic_variables()) {
    //cout << "    " << b << endl;
    //assert(!elem(b, tab.basic_variables));
  //}

  //sanity_check(tab);

  //for (int r = 0; r < tab.num_rows(); r++) {
    //for (int c = 0; c < tab.num_cols(); c++) {
      //cout << (tab.get_entry(r, c)) << " ";
    //}
    //cout << endl;
  //}

  //while (can_improve(tab)) {
    //int next_pivot_col = pick_pivot_col(tab);
    //cout << "next pivot col = " << next_pivot_col << endl;
    //bool unbounded = true;
    //for (int r = 0; r < tab.num_rows(); r++) {
      //if (tab.get_entry(r, next_pivot_col) > 0) {
        //unbounded = false;
        //break;
      //}
    //}

    //if (unbounded) {
      //return {0, LP_RESULT_TYPE_UNBOUNDED};
    //}

    //int pivot_row = pick_pivot_row(next_pivot_col, tab);
    //cout << "Pivot row = " << pivot_row << endl;

    //assert(pivot_row >= 0);

    //int next_pivot_row = pivot_row;
    //cout << "--- Next non basic variable to convert: " << next_pivot_col << endl;

    //auto pivot_val = tab.get_entry(next_pivot_row, next_pivot_col);
    //cout << "pivot val = " << pivot_val << endl;

    //assert(pivot_val != 0);

    //tab.maximum = tab.maximum - ((tab.objective_coeff(next_pivot_col)) * ((tab.const_coeff(next_pivot_row)) / pivot_val));
    //cout << "Pivot val = " << pivot_val << endl;
    //cout << "New max   = " << (tab.maximum) << endl;
    //tab.scale_row((value(1) / pivot_val), pivot_row);

    //for (int r = 0; r < tab.num_rows(); r++) {
      //if (r != next_pivot_row) {
        //vector<value> muls;
        //for (int c = 0; c < tab.num_cols(); c++) {
          //muls.push_back(tab.get_entry(r, c) * tab.get_entry(pivot_row, c));
        //}
        //tab.subtract_row(muls, r);
      //}
    //}

    //tab.exchange(next_pivot_col, next_pivot_row);
    //tab.print(cout);
  //}

  //cout << "#### Final tableau" << endl;
  //tab.print(cout);
  //for (int r = 1; r < tab.num_rows(); r++) {
    //bool contains_one = false;
    //for (int c = 0; c < tab.num_cols(); c++) {
      //if (tab.get_entry(r, c) == 0) {
        //contains_one = true;
        //break;
      //}
    //}
    //if (contains_one) {
      //basis.insert(r);
    //}
  //}


  //return {tab.maximum, LP_RESULT_TYPE_OPTIMAL};
}

linear_constraint eq(const linear_expr& e) {
  return {LINEAR_CONSTRAINT_TYPE_EQ, e};
}

linear_constraint geq(const linear_expr& e) {
  return {LINEAR_CONSTRAINT_TYPE_GEQ, e};
}

linear_constraint leq(const linear_expr& e) {
  return {LINEAR_CONSTRAINT_TYPE_LEQ, e};
}

value evaluate(const point& p, const linear_expr& lc) {
  assert(p.dimension() == lc.dimension());
  value v = 0;
  for (int i = 0; i < p.dimension(); i++) {
    v = v + p[i] * lc[i];
  }
  v = v + lc.constant();
  return v;
}

bool sat(const point& p, const linear_constraint& lc) {
  //assert(p.dimension() == lc.dimension());
  //value v = 0;
  //for (int i = 0; i < p.dimension(); i++) {
    //v += p[i] + lc[i];
  //}
  //v += lc.constant();
  return evaluate(p, lc.expr) >= 0;
}

bool tight(const point& p, const linear_constraint& lc) {
  assert(p.dimension() == lc.dimension());
  auto r = evaluate(p, lc.expr);
  cout << "r = " << r << endl;
  return r == 0;
}

void constraint_set_test() {
  linear_expr sum(1);
  sum.set_coeff(0, value(1));

  linear_expr lc(1);
  lc.set_coeff(0, value(-1));
  lc.set_const(value(5));

  linear_constraint constraint = geq(lc);

  cout << "C = " << constraint << endl;

  assert(sat({0}, constraint));
  assert(tight({5}, constraint));
  assert(!tight({3}, constraint));
}

void basic_test() {
  linear_expr sum(1);
  sum.set_coeff(0, value(1));

  linear_expr lc(1);
  lc.set_coeff(0, value(-1));
  lc.set_const(value(5));

  vector<linear_constraint> constraints{geq(lc)};
  lp_result result = maximize(sum, constraints);

  assert(result.val == 5);
  assert(result.tp == LP_RESULT_TYPE_OPTIMAL);

  cout << "BASIC TEST Passed" << endl;
  assert(false);
}

//void ft_test() {
  //linear_expr sum(1);
  //sum.set_coeff(0, value(1));

  //linear_expr lc(1);
  //lc.set_coeff(0, value(-2));
  //lc.set_const(value(60));

  //vector<linear_expr> constraints{lc};
  //lp_result result = maximize(sum, constraints);

  //assert(result.val == 30);
  //assert(result.tp == LP_RESULT_TYPE_OPTIMAL);

  //cout << "Non-unity coefficient test passed" << endl;
//}

//void unbounded_test() {
  //linear_expr sum(1);
  //sum.set_coeff(0, value(1));

  //linear_expr lc(1);
  //lc.set_coeff(0, value(1));
  //lc.set_const(value(5));

  //vector<linear_expr> constraints{lc};
  //lp_result result = maximize(sum, constraints);

  //assert(result.tp == LP_RESULT_TYPE_UNBOUNDED);

  //cout << "Unbounded test passed" << endl;
//}

//void no_solution_test() {
  //cout << "------ Starting no_solution_test" << endl;
  //linear_expr sum(1);
  //sum.set_coeff(0, value(1));

  //linear_expr lc(1);
  //lc.set_coeff(0, value(1));
  //lc.set_const(value(5));

  //linear_expr cc(1);
  //cc.set_coeff(0, value(-1));
  //cc.set_const(value(-6));

  //cout << "cc=  " << cc << endl;
  //cout << "lc=  " << lc << endl;

  //lp_result result = maximize(sum, {lc, cc});
  //cout << "result max = " << result.val << endl;
  //assert(result.tp == LP_RESULT_TYPE_INFEASIBLE);

  //cout << "No solution test passed" << endl;
//}

void phase_1_test() {
  tableau tab(6, 11);
  tab(0, 0) = 1;
  tab(0, 8) = 1;
  tab(0, 9) = 1;

  tab(1, 1) = 1;
  tab(1, 2) = -5;
  tab(1, 3) = 1;
  tab(1, 4) = 1;

  tab(2, 2) = 3;
  tab(2, 3) = -1;
  tab(2, 4) = -1;
  tab(2, 5) = 1;
  tab(2, 9) = -1;
  tab(2, 9) = -1;

  tab(3, 2) = 1;
  tab(3, 3) = 2;
  tab(3, 4) = -1;
  tab(3, 6) = 1;
  tab(3, 9) = -1;
  tab(3, 9) = -2;

  tab(4, 2) = 2;
  tab(4, 3) = 1;
  tab(4, 7) = 1;
  tab(4, 9) = 2;

  tab(5, 2) = 1;
  tab(5, 3) = 1;
  tab(5, 8) = 1;
  tab(5, 9) = 1;

  cout << tab << endl;
  assert(false);
}

int main() {
  constraint_set_test();
  phase_1_test();
  basic_test();
  //ft_test();
  //unbounded_test();
  //no_solution_test();
}

