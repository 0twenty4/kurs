#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <iostream>
#include <algorithm>

struct NODE {
	double x, y;
	int global_num = -1;
	
	int formula_num = -1;

	NODE(double x, double y, int global_num) :
		x(x), y(y), global_num(global_num) {}

	NODE(double x, double y) :
		x(x), y(y) {}

	NODE() {}

	bool operator<(const NODE& node) const {
		return global_num < node.global_num;
	}
};

struct ELEM {
	std::vector<NODE> nodes;

	NODE mass_center;

	int formula_num = -1;

	std::vector<std::vector<double>> local_mat;
	std::vector<double> local_vec;

	ELEM(std::vector<NODE> nodes, int formula_num, int n) :
		local_mat(n, std::vector<double>(n, 0)), local_vec(n, 0),
		nodes(nodes), formula_num(formula_num) {}

	ELEM(std::vector<NODE> nodes) :
		nodes(nodes) {}
};

struct AREA {
	std::vector<ELEM> grid;
	std::vector<NODE> nodes;

	std::vector<ELEM> first_bound_cond, second_bound_cond, third_bound_cond;
};

struct MAT {
	int n;
	std::vector<int> ig, jg;
	std::vector<double> ggl, di;

	MAT(int n) : n(n) {}
};

struct SLAE {
	MAT mat;
	std::vector<double> vec;
	std::vector<double> sol;
	double epsilon;

	SLAE(int n, double epsilon) : mat(n), vec(n), sol(n), epsilon(epsilon) {}

	SLAE(MAT mat,int n) :
		mat(mat), vec(n), sol(n) {}

	void clear_sol() {
		for (int i = 0; i < mat.n; i++)
			sol[i] = 0;
	}
};

void input(AREA& area, std::string nodes_file, std::string elems_file,
	std::string first_bound_cond, std::string second_bound_cond_file,
	std::string third_bound_cond_file) {

	std::ifstream in(nodes_file);

	double x, y;
	for (int global_num = 0; in >> x >> y; global_num++) {
		NODE node(x, y, global_num);
		area.nodes.push_back(node);
	}

	in.close();
	in.open(elems_file);

	int node1_global_num, node2_global_num, node3_global_num, formula_num;

	while (in >> node1_global_num >> node2_global_num >> node3_global_num >> formula_num) {

		ELEM elem({ area.nodes[node1_global_num],area.nodes[node2_global_num],area.nodes[node3_global_num] },
			formula_num, 3);

		std::sort(elem.nodes.begin(), elem.nodes.end());
		
		elem.mass_center.x = (elem.nodes[0].x + elem.nodes[1].x + elem.nodes[2].x) / 3;
		elem.mass_center.y = (elem.nodes[0].y + elem.nodes[1].y + elem.nodes[2].y) / 3;

		area.grid.push_back(elem);
	}
	in.close();

	in.open(first_bound_cond);


	while (in >> node1_global_num >> node2_global_num >> formula_num) {

		ELEM edge({ area.nodes[node1_global_num], area.nodes[node2_global_num] }, formula_num, 0);
		area.first_bound_cond.push_back(edge);
	}

	in.close();

	in.open(second_bound_cond_file);

	while (in >> node1_global_num >> node2_global_num >> formula_num) {
		ELEM edge({ area.nodes[node1_global_num], area.nodes[node2_global_num] }, formula_num, 2);
		area.second_bound_cond.push_back(edge);
	}


	in.close();
	in.open(third_bound_cond_file);

	while (in >> node1_global_num >> node2_global_num >> formula_num) {
		ELEM edge({ area.nodes[node1_global_num], area.nodes[node2_global_num] }, formula_num, 2);
		area.third_bound_cond.push_back(edge);
	}

	in.close();
}

double u(NODE node) {
	return node.x + 6 * node.y - 2;
}

double calc_detD(ELEM& elem) {
	return (elem.nodes[1].x - elem.nodes[0].x) * (elem.nodes[2].y - elem.nodes[0].y) -
		(elem.nodes[2].x - elem.nodes[0].x) * (elem.nodes[1].y - elem.nodes[0].y);
}

double calc_mes_edge(ELEM& edge) {

	return sqrt(pow(edge.nodes[1].x - edge.nodes[0].x, 2) + pow(edge.nodes[1].y - edge.nodes[0].y, 2));

}

double mass_center_num_val(std::vector<double>& slae_sol, ELEM& elem) {
	double elem_detD = calc_detD(elem);

	ELEM S23({ elem.nodes[1],elem.nodes[2],elem.mass_center });
	ELEM S31({ elem.nodes[2], elem.nodes[0],elem.mass_center });
	ELEM S12({ elem.nodes[0],elem.nodes[1],elem.mass_center });

	double L1 = calc_detD(S23) / elem_detD;
	double L2 = calc_detD(S31) / elem_detD;
	double L3 = calc_detD(S12) / elem_detD;

	return slae_sol[elem.nodes[0].global_num] * L1 +
		slae_sol[elem.nodes[1].global_num] * L2 +
		slae_sol[elem.nodes[2].global_num] * L3;
}

void output(AREA& area, std::vector<double> slae_sol, std::string res_file) {
	std::ofstream out(res_file);

	for (int i = 0; i < slae_sol.size(); i++)
		out << slae_sol[i] << " " << u(area.nodes[i]) << std::endl;

	out << std::endl;

	for (int i = 0; i < area.grid.size(); i++)
		out << mass_center_num_val(slae_sol, area.grid[i]) << " " << u(area.grid[i].mass_center) << std::endl;
}

void portrait(MAT & mat, AREA area) {
	std::vector< std::set<int>> list(area.nodes.size());

	for (int elem = 0; elem < area.grid.size(); elem++) {
		for (int node1 = 2; node1 >= 0; node1--) {
			for (int node2 = node1-1; node2 >= 0; node2--) {
				list[area.grid[elem].nodes[node1].global_num].insert(area.grid[elem].nodes[node2].global_num);
			}
		}
	}

	mat.ig.resize(area.nodes.size() + 1);
	mat.ig[0] = 0;

	for (int i = 1; i < area.nodes.size() + 1; i++) {
		mat.ig[i] = mat.ig[i - 1] + list[i - 1].size();
	}

	mat.ggl.resize(mat.ig.back());
	mat.jg.resize(mat.ggl.size());
	mat.di.resize(area.nodes.size());

	for (int node = 0, k = 0; node < list.size(); node++) {
		for (auto conn_node = list[node].begin(); conn_node != list[node].end(); ++conn_node, k++)
			mat.jg[k] = *conn_node;
	}
}

double gamma(int formula_num) {
	switch (formula_num) {
	case(1):
		return 5;
		break;
	case(2):
		return 0;
		break;
	}
}

double lambda(int formula_num, NODE node) {
	switch (formula_num) {
	case(1):
		return node.y;
		break;
	case(2):
		return node.y;
		break;
	}
		
}

double f(int formula_num, NODE node) {
	switch (formula_num) {
	case(1):
		return 5*node.x+30*node.y-16;
		break;
	case(2):
		return -6;
		break;
	}
}

double ug(int formula_num, NODE node) {
	switch (formula_num) {
	case(1):
		return 6*node.y+2;
		break;
	}
}

double beta(int formula_num) {
	switch (formula_num) {
	case(1):
		return 1;
		break;
	}
}

double ubeta(int formula_num, NODE node) {
	switch (formula_num) {
	case(1):
		return 7*node.y+2;
		break;
	}
}

double teta(int formula_num, NODE node) {
	switch (formula_num) {
	case(1):
		return -6*node.y;
		break;
	case(2):
		return -node.y;
		break;
	case(3):
		return 6*node.y;
		break;
	}
}

void solve_local(ELEM & elem) {
	double detD = calc_detD(elem);
	std::vector<std::vector<double>> coefs =

	{ {(elem.nodes[1].x * elem.nodes[2].y - elem.nodes[2].x * elem.nodes[1].y) / detD,
		(elem.nodes[1].y - elem.nodes[2].y) / detD,
		(elem.nodes[2].x - elem.nodes[1].x) / detD},

		{(elem.nodes[2].x * elem.nodes[0].y - elem.nodes[0].x * elem.nodes[2].y) / detD,
		(elem.nodes[2].y - elem.nodes[0].y) / detD,
		(elem.nodes[0].x - elem.nodes[2].x) / detD},

		{(elem.nodes[0].x * elem.nodes[1].y - elem.nodes[1].x * elem.nodes[0].y) / detD,
		(elem.nodes[0].y - elem.nodes[1].y) / detD,
		(elem.nodes[1].x - elem.nodes[0].x) / detD}

	};

	NODE mid_node_1(
		(elem.nodes[0].x + elem.nodes[1].x) / 2.,
		(elem.nodes[0].y + elem.nodes[1].y) / 2.
	);

	NODE mid_node_2{
		(elem.nodes[1].x + elem.nodes[2].x) / 2.,
		(elem.nodes[1].y + elem.nodes[2].y) / 2.,
	};

	NODE mid_node_3{
		(elem.nodes[0].x + elem.nodes[2].x) / 2.,
		(elem.nodes[0].y + elem.nodes[2].y) / 2.,
	};

	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			elem.local_mat[i][j] =

				(i == j ? gamma(elem.formula_num) * (abs(detD) / 12) :
					gamma(elem.formula_num) * (abs(detD) / 24)) +

				(abs(detD) / 6) * (coefs[i][1] * coefs[j][1] + coefs[i][2] * coefs[j][2]) *
				(lambda(elem.formula_num, mid_node_1)
					+ lambda(elem.formula_num, mid_node_2)
					+ lambda(elem.formula_num, mid_node_3));
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int k = 0; k < 3; k++) {
			elem.local_vec[i] += i == k ?
				f(elem.formula_num, elem.nodes[k]) * (abs(detD) / 12) :
				f(elem.formula_num, elem.nodes[k]) * (abs(detD) / 24);
		}
	}
}

int find_ind_dichotomy(std::vector<int>&arr, int start, int end, int goal) {
	int i;

	while (start<=end) {

		i = start + (end - start) / 2.;

		if (arr[i] == goal)
			return i;

		if (arr[i] < goal)
			start = i + 1;
		else
			end = i-1;
	}

	return -1;
}

void add_to_global(SLAE& slae, ELEM & elem) {

	for (int i = 0; i < elem.local_mat.size(); i++) {
		slae.mat.di[elem.nodes[i].global_num] += elem.local_mat[i][i];
		slae.vec[elem.nodes[i].global_num] += elem.local_vec[i];

		for (int j = 0; j < i; j++) {
			int ggl_ij = find_ind_dichotomy(slae.mat.jg,
				slae.mat.ig[elem.nodes[i].global_num],
				slae.mat.ig[elem.nodes[i].global_num + 1]-1,
				elem.nodes[j].global_num);

				slae.mat.ggl[ggl_ij] += elem.local_mat[i][j];

		}
	}
}

void first_bound_cond(SLAE& slae, ELEM& edge) {
	for (int node = 0; node < 2;  node++) {
		slae.mat.di[edge.nodes[node].global_num] = 1;
		slae.vec[edge.nodes[node].global_num] = ug(edge.formula_num, edge.nodes[node]);

		for (int k = 0; k < slae.mat.ig[edge.nodes[node].global_num + 1] - slae.mat.ig[edge.nodes[node].global_num]; k++) {
			slae.vec[slae.mat.jg[slae.mat.ig[edge.nodes[node].global_num]+k]] +=
				slae.vec[edge.nodes[node].global_num] * (-slae.mat.ggl[slae.mat.ig[edge.nodes[node].global_num] + k]);

			slae.mat.ggl[slae.mat.ig[edge.nodes[node].global_num] + k] = 0;
		}

		for (int i = edge.nodes[node].global_num + 1; i < slae.mat.n; i++) {
			int ggl_inodegm = find_ind_dichotomy(slae.mat.jg,
				slae.mat.ig[i],
				slae.mat.ig[i + 1] - 1,
				edge.nodes[node].global_num);

			if (ggl_inodegm != -1) {
				slae.vec[i] += slae.vec[edge.nodes[node].global_num] * (-slae.mat.ggl[ggl_inodegm]);
				slae.mat.ggl[ggl_inodegm] = 0;
			}
		}
	}
}

void second_bound_cond(SLAE& slae, ELEM& edge) {
	std::vector<std::vector<double>> coefs = {
		{2.,1.},
		{1.,2.}
	};

	double mes_edge = calc_mes_edge(edge);

	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < 2; k++) {
			edge.local_vec[i] += (mes_edge / 6.) * coefs[i][k] * teta(edge.formula_num, edge.nodes[k]);
		}
	}

	add_to_global(slae, edge);
}

void third_bound_cond(SLAE& slae, ELEM& edge) {
	double mes_edge = calc_mes_edge(edge);

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			edge.local_mat[i][j] = i == j ? 
				beta(edge.formula_num) * mes_edge / 3. :
				beta(edge.formula_num) * mes_edge / 6.;
		}
	}

	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < 2; k++) {
			edge.local_vec[i] += edge.local_mat[i][k] * ubeta(edge.formula_num, edge.nodes[k]);
		}
	}

	add_to_global(slae, edge);
}

MAT factor(MAT& mat) {
	MAT res = mat;
	double sum;
	int found_num1, found_num2;
	for (int i = 0; i < res.n; i++) {
		sum = 0;
		found_num1 = found_num2 = 0;
		for (int j = 0; j < i; j++) {
			int ggl_ij = find_ind_dichotomy(res.jg,
				res.ig[i]+found_num1, res.ig[i + 1] - 1, j);

			if (ggl_ij != -1) {
				for (int k = 0; k < j; k++) {
					int ggl_ik = find_ind_dichotomy(res.jg,
						res.ig[i]+found_num2, res.ig[i + 1] - 1, k);
					int ggl_jk = find_ind_dichotomy(res.jg,
						res.ig[j]+found_num2, res.ig[j + 1] - 1, k);

					if (ggl_ik != -1 and ggl_jk != -1)
						sum += res.ggl[ggl_ik] * res.ggl[ggl_jk];
					
					found_num2++;
				}

				res.ggl[ggl_ij] = (1 / res.di[j]) * (res.ggl[ggl_ij] - sum);
				found_num1++;
			}
			sum = 0;
		}

		found_num1 = 0;
		for (int k = 0; k < i; k++) {
			int ggl_ik = find_ind_dichotomy(res.jg,
				res.ig[i]+found_num1, res.ig[i + 1] - 1, k);

			if (ggl_ik != -1) {
				sum += res.ggl[ggl_ik] * res.ggl[ggl_ik];
				found_num1++;
			}
		}

		res.di[i] = sqrt(res.di[i] - sum);

	}

	return res;
}

void solve_L(SLAE& slae) {
	int found_num;
	for (int i = 0; i < slae.mat.n; i++) {
		slae.sol[i] = slae.vec[i];

		found_num = 0;
		for (int j = 0; j < i; j++) {
			int ggl_ij = find_ind_dichotomy(slae.mat.jg,
				slae.mat.ig[i]+ found_num, slae.mat.ig[i + 1] - 1,
				j);

			if (ggl_ij != -1) {
				slae.sol[i] -= slae.mat.ggl[ggl_ij] * slae.sol[j];
				found_num++;
			}
		}
		slae.sol[i] /=slae.mat.di[i];
	}
}

void solve_Lt(SLAE& slae) {
	int found_num;
	for (int j = slae.mat.n - 1; j >= 0; j--) {
		slae.sol[j] /= slae.mat.di[j];

		found_num = 0;
		for (int i = j - 1; i >= 0; i--) {
			int ggl_ji = find_ind_dichotomy(slae.mat.jg,
				slae.mat.ig[j], slae.mat.ig[j + 1] - 1-found_num,
				i);

			if (ggl_ji != -1) {
				slae.sol[i] -= slae.mat.ggl[ggl_ji] * slae.sol[j];
				found_num++;
			}
		}
	}
}

void solve_SLAE_direct(SLAE& slae) {
	solve_L(slae);
	solve_Lt(slae);
}

std::vector<double> mult_mat_CSR_vec(MAT mat, std::vector<double> vec) {
	std::vector<double> res(mat.n);
	int found_num;
	for (int i = 0; i < mat.n; i++) {
		found_num = 0;
		res[i] += mat.di[i] * vec[i];
		for (int j = 0; j < i; j++) {
			int ggl_ij = find_ind_dichotomy(mat.jg,
				mat.ig[i]+found_num, mat.ig[i + 1] - 1, j);

			if (ggl_ij != -1) {
				res[i] += mat.ggl[ggl_ij] * vec[j];
				res[j] += mat.ggl[ggl_ij] * vec[i];
				found_num++;
			}
		}
	}
	return res;
}

std::vector<double> vec_sum(std::vector<double> vec1, std::vector<double> vec2, int sign) {
	std::vector<double> res(vec1.size());

	for (int i = 0; i < res.size(); i++) {
		res[i] = sign == 1 ? vec1[i] + vec2[i] : vec1[i] - vec2[i];
	}

	return res;
}

double scal_mult(std::vector<double> vec1, std::vector<double> vec2) {
	double res = 0;

	for (int i = 0; i < vec1.size(); i++)
		res += vec1[i] * vec2[i];

	return res;
}

std::vector<double> mult_scal_vec(double scal, std::vector<double> vec) {
	std::vector<double> res(vec.size());

	for (int i = 0; i < vec.size(); i++)
		res[i] = scal * vec[i];

	return res;

}

void solve_SLAE_CGM(SLAE& slae, int max_iter) {
	std::vector<double> r_prev = vec_sum(slae.vec, mult_mat_CSR_vec(slae.mat, slae.sol),-1);
	std::vector<double> r_cur = r_prev;
	SLAE Mr_prev_slae(factor(slae.mat), slae.mat.n);
	SLAE Mr_cur_slae = Mr_prev_slae;

	Mr_prev_slae.vec = r_prev;

	solve_SLAE_direct(Mr_prev_slae);

	std::vector<double> z = Mr_prev_slae.sol;

	std::vector<double> Az;
	double a, b;

	double f_norm = sqrt(scal_mult(slae.vec, slae.vec));
	double r_cur_norm;

	int iter = 0;
	while(true) {
		iter++;

		Az = mult_mat_CSR_vec(slae.mat, z);

		a = scal_mult(Mr_prev_slae.sol, r_prev) / scal_mult(Az, z);

		slae.sol = vec_sum(slae.sol, mult_scal_vec(a, z), 1);

		r_cur = vec_sum(r_prev, mult_scal_vec(a, Az), -1);

		r_cur_norm = sqrt(scal_mult(r_cur, r_cur));

		std::cout << r_cur_norm / f_norm << std::endl;

		if (r_cur_norm / f_norm < slae.epsilon or iter == max_iter)
			break;

		Mr_cur_slae.vec = r_cur;
		Mr_cur_slae.clear_sol();
		solve_SLAE_direct(Mr_cur_slae);

		b = scal_mult(Mr_cur_slae.sol, r_cur) / scal_mult(Mr_prev_slae.sol, r_prev);

		z = vec_sum(Mr_cur_slae.sol, mult_scal_vec(b, z), 1);
		
		r_prev = r_cur;

		Mr_prev_slae.sol = Mr_cur_slae.sol;
	} 
}

int main() {
	AREA area;

	input(area, "nodes.txt", "elems.txt", "first_bound_cond.txt", "second_bound_cond.txt", "third_bound_cond.txt");

	MAT mat(area.nodes.size());

	portrait(mat, area);

	SLAE slae(mat, mat.n);

	slae.epsilon = 1e-12;

	for (ELEM& elem : area.grid) {
		solve_local(elem);
		add_to_global(slae, elem);
	}

	for (ELEM& edge : area.third_bound_cond)
		third_bound_cond(slae, edge);

	for (ELEM& edge : area.second_bound_cond)
		second_bound_cond(slae, edge);

	for (ELEM& edge : area.first_bound_cond)
		first_bound_cond(slae, edge);

	for (auto& elem : area.grid) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				std::cout << elem.local_mat[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "\n\n";

	}

	for (auto& elem:area.third_bound_cond) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				std::cout << elem.local_mat[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "\n\n";

	}

	for (auto& di : slae.mat.di)
		std::cout << di << " ";

	std::cout << std::endl;
	
	for (auto& ggl : slae.mat.ggl)
		std::cout << ggl << " ";

	std::cout << std::endl;

	for (auto& ig : slae.mat.ig)
		std::cout << ig << " ";

	std::cout << std::endl;

	for (auto& jg : slae.mat.jg)
		std::cout << jg << " ";

	std::cout << std::endl;

	for (auto& vec : slae.vec)
		std::cout << vec << " ";

	std::cout << std::endl;

	solve_SLAE_CGM(slae, 1000000);

	std::vector<double> mass_center_num_vals(area.grid.size());

	for (int i=0; i<area.grid.size(); i++)
		mass_center_num_vals[i] = mass_center_num_val(slae.sol, area.grid[i]);

	output(area, slae.sol, "res.txt");

}
