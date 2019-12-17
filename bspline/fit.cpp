#include "fit.hpp"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <algorithm>
#include <glm/glm.hpp>
namespace {
	constexpr float eps = 1e-7f;
}

std::vector<vec3> line_interplation(
	const std::vector<vec3>& data,
	int p,
	int m,
	const std::vector<float>& ts,
	const Bspline& b) {
	int n = data.size();
	Eigen::MatrixXf M(n, m + 1);
	Eigen::MatrixXf D(n, 3);
	for (int i = 0; i < n; i++) {
		auto nik = b.Nik(p + 1, ts[i]);
		for (int j = 0; j < m+1; j++)
			M(i, j) = nik[j];
		for (int j = 0; j < 3; j++)
			D(i, j) = data[i][j];
	}
	Eigen::MatrixXf P = M.lu().solve(D);
	std::vector<vec3> res(m + 1);
	for (int i = 0; i < m+1; i++) {
		res[i].x = P(i, 0);
		res[i].y = P(i, 1);
		res[i].z = P(i, 2);
	}
	return res;
}

bool Fitter::interplation(const std::vector<std::vector<vec3>>& data, int p,
	int q, int m, int n,
	std::vector<std::vector<vec3>>* ps,
	std::vector<float>* us, std::vector<float>* vs,
	Bspline_Surface** bs) {
	*bs = get_parameter(data, p, q, m, n, us, vs);
	int num_col = vs->size();
	int num_row = us->size();
	std::vector<std::vector<vec3>> Q(num_col, std::vector<vec3>(m + 1));
	for (int col = 0; col < num_col; col++) {
		std::vector<vec3> col_data(num_row);
		for (int row = 0; row < num_row; row++) {
			col_data[row] = data[row][col];
		}
		Q[col] = line_interplation(col_data, p, m, *us, (*bs)->ubs);
	}
	for (int row = 0; row < m + 1; row++) {
		std::vector<vec3> row_data(num_col);
		for (int col = 0; col < num_col; col++) {
			row_data[col] = Q[col][row];
		}
		(*ps)[row] = line_interplation(row_data, q, n, *vs, (*bs)->vbs);
	}
	return true;
	/*int un = us->size();
	int vn = vs->size();
	std::vector<std::vector<vec3>> Q(un, std::vector<vec3>(vn));
	for (int col = 0; col < vn; col++) {
		Eigen::MatrixXf M(un, m + 1);
		Eigen::MatrixXf D(un, 3);
		for (int i = 0; i < un; i++) {
			auto nik = (*bs)->ubs.Nik(p + 1, (*us)[i]);
			for (int j = 0; j < nik.size(); j++)
				M(i, j) = nik[j];
			for (int j = 0; j < 3; j++)
				D(i, j) = data[i][col][j];
		}
		Eigen::MatrixXf res = M.lu().solve(D);
		for (int i = 0; i < un; i++) {
			Q[i][col].x = res(i, 0);
			Q[i][col].y = res(i, 1);
			Q[i][col].z = res(i, 2);
		}
	}
	printf("\nD\n");
	for (int row = 0; row < un; row++) {
		Eigen::MatrixXf M(vn, n + 1);
		Eigen::MatrixXf D(vn, 3);
		for (int i = 0; i < vn; i++) {
			auto nik = (*bs)->vbs.Nik(q + 1, (*vs)[i]);
			for (int j = 0; j < nik.size(); j++)
				M(i, j) = nik[j];
			for (int j = 0; j < 3; j++)
				D(i, j) = Q[row][i][j];
		}
		Eigen::MatrixXf res = M.lu().solve(D);
		for (int i = 0; i < vn; i++) {
			(*ps)[row][i].x = res(i, 0);
			(*ps)[row][i].y = res(i, 1);
			(*ps)[row][i].z = res(i, 2);
		}
		auto d = M * res;
		for (int i = 0; i < vn; i++) {
			printf("%f %f %f  ", d(i, 0), d(i, 1), d(i, 2));
		}
		printf("\n");
	}
	printf("\nQ\n");
	for (auto&& r : Q) {
		for (auto&& c : r) {
			printf("%f %f %f  ", c.x, c.y, c.z);
		}
		printf("\n");
	}
	printf("\nR\n");
	for (int i = 0; i < un; i++) {
		for (int j = 0; j < vn; j++) {
			auto v = (*bs)->vbs((*ps)[i], (*vs)[j]);
			printf("%f %f %f  ", v.x, v.y, v.z);
		}
		printf("\n");
	}*/
	return true;
}

std::vector<vec3> line_approximation(
	const std::vector<vec3>& data,
	int p,
	int h,
	const std::vector<float>& ts,
	const Bspline& b) {
	int n = ts.size()-1;
	Eigen::MatrixXf N(n - 1, h - 1);
	Eigen::MatrixXf Q(h - 1, 3);
	Eigen::MatrixXf Qk(n - 1, 3);
	vec3 D0 = data.front();
	vec3 Dn = data.back();
	for (int i = 1; i <= n - 1; i++) {
		auto&& nik = b.Nik(p + 1, ts[i]);
		auto qk = data[i] - nik.front() * D0 - nik.back() * Dn;
		for (int j = 0; j < 3; j++) Qk(i - 1, j) = qk[j];
		for (int j = 1; j <= h - 1; j++) {
			N(i - 1, j - 1) = nik[j];
		}
	}
	Q = (Qk.transpose() * N).transpose();
	Eigen::MatrixXf P = (N.transpose() * N).lu().solve(Q);
	std::vector<vec3> res(h + 1);
	res[0] = D0;
	res[h] = Dn;
	for (int i = 1; i <= h - 1; i++) {
		res[i].x = P(i - 1,0);
		res[i].y = P(i - 1,1);
		res[i].z = P(i - 1,2);
	}
	return res;
}

bool Fitter::approximation(const std::vector<std::vector<vec3>>& data, int p,
	int q, int m, int n,
	std::vector<std::vector<vec3>>* ps,
	std::vector<float>* us, std::vector<float>* vs,
	Bspline_Surface** bs) {
	*bs = get_parameter(data, p, q, m, n, us, vs);
	int num_col = vs->size();
	int num_row = us->size();
	std::vector<std::vector<vec3>> tmp_data(num_col, std::vector<vec3>(m + 1));
	for (int col = 0; col < vs->size(); col++) {
		std::vector<vec3> col_data(num_row);
		for (int row = 0; row < num_row; row++) {
			col_data[row] = data[row][col];
		}
		tmp_data[col] = line_approximation(col_data, p, m, *us, (*bs)->ubs);
	}
	for (int row = 0; row < m + 1; row++) {
		std::vector<vec3> row_data(num_col);
		for (int col = 0; col < num_col; col++) {
			row_data[col] = tmp_data[col][row];
		}
		(*ps)[row] = line_approximation(row_data, q, n, *vs, (*bs)->vbs);
	}
	return true;
}

Bspline_Surface*
Fitter::get_parameter(const std::vector<std::vector<vec3>>& data, int p, int q,
	int m, int n, std::vector<float>* us,
	std::vector<float>* vs) const {
	int h = data.size();
	int w = data[0].size();
	std::vector<std::vector<float>> uvs(h, std::vector<float>(w, 0));
	std::vector<float> knots_u, knots_v;
	auto f_uniform = [](const std::vector<vec3>& points) {
		int n = points.size();
		std::vector<float> ts(n);
		float delta = 1.f / (n - 1);
		float sum = 0;
		for (int i = 0; i < n; i++) {
			ts[i] = sum;
			sum += delta;
		}
		ts[n - 1] = 1.f;
		return ts;
	};
	auto f_chordlength = [](const std::vector<vec3>& points) {
		int n = points.size();
		std::vector<float> arcs;
		float L = 0;
		for (int i = 1; i < n; i++) {
			float Li = glm::length(points[i] - points[i - 1]);
			arcs.emplace_back(Li);
			L += Li;
		}
		std::vector<float> ts(n, 0);
		float sum = 0;
		for (int i = 1; i < n; i++) {
			sum += arcs[i - 1];
			ts[i] = sum / L;
		}
		ts[n - 1] = 1.f;
		return ts;
	};
	auto f_centripetal = [](const std::vector<vec3>& points, float alpha) {
		int n = points.size();
		std::vector<float> arcs;
		float L = 0;
		for (int i = 1; i < n; i++) {
			float Li = std::pow(glm::length(points[i] - points[i - 1]), alpha);
			arcs.emplace_back(Li);
			L += Li;
		}
		std::vector<float> ts(n, 0);
		float sum = 0;
		for (int i = 1; i < n; i++) {
			sum += arcs[i - 1];
			ts[i] = sum / L;
		}
		ts[n - 1] = 1.f;
		return ts;
	};
	auto f_uniform_knots = [](const std::vector<float>& us, int p, int n) {
		int m = p + n + 1;
		std::vector<float> ts(m + 1);
		for (int i = 0; i <= p; i++)
			ts[i] = 0;
		for (int i = m - p; i <= m; i++)
			ts[i] = 1.f + eps;
		float delta = 1.f / (n - p + 1);
		float sum = delta;
		for (int i = p + 1; i <= n; i++) {
			ts[i] = sum;
			sum += delta;
		}
		return ts;
	};
	auto f_average_knots = [](const std::vector<float>& us, int p, int n) {
		int m = p + n + 1;
		std::vector<float> ts(m + 1);
		for (int i = 0; i <= p; i++)
			ts[i] = 0;
		for (int i = m - p; i <= m; i++)
			ts[i] = 1.f + eps;
		int step = us.size() - 1 - n + p;
		float inv_p = 1.f / step;
		for (int j = 1; j <= n - p; j++) {
			float sum = 0.f;
			for (int i = j; i <= j + step - 1; i++) {
				sum += us[i];
			}
			ts[j + p] = inv_p * sum;
		}
		return ts;
	};
	if (parameter_method != universal) {
		std::vector<vec3> points(h);
		for (int j = 0; j < w; j++) {
			for (int i = 0; i < h; i++) {
				points[i] = data[i][j];
			}
			std::vector<float> ts;
			if (parameter_method == uniformly_space) {
				ts = f_uniform(points);
			}
			else if (parameter_method == chordlength) {
				ts = f_chordlength(points);
			}
			else if (parameter_method == centripetal) {
				ts = f_centripetal(points, alpha);
			}
			else {
				assert(0);
			}
			for (int i = 0; i < h; i++) {
				uvs[i][j] = ts[i];
			}
		}
		for (int i = 0; i < h; i++) {
			float sum = 0;
			for (int j = 0; j < w; j++)
				sum += uvs[i][j];
			(*us)[i] = sum / w;
		}

		points.resize(w);
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				points[j] = data[i][j];
			}
			std::vector<float> ts;
			if (parameter_method == uniformly_space) {
				ts = f_uniform(points);
			}
			else if (parameter_method == chordlength) {
				ts = f_chordlength(points);
			}
			else if (parameter_method == centripetal) {
				ts = f_centripetal(points, alpha);
			}
			else {
				assert(0);
			}
			for (int j = 0; j < w; j++) {
				uvs[i][j] = ts[j];
			}
		}
		for (int j = 0; j < w; j++) {
			float sum = 0;
			for (int i = 0; i < h; i++)
				sum += uvs[i][j];
			(*vs)[j] = sum / h;
		}

		if (knot_generation == uniform) {
			knots_u = f_uniform_knots(*us, p, m);
			knots_v = f_uniform_knots(*vs, q, n);
		}
		else if (knot_generation == average) {
			knots_u = f_average_knots(*us, p, m);
			knots_v = f_average_knots(*vs, q, n);
		}
		else {
			assert(0);
		}
	}
	else {
		/*
		universal method can not be used to approximate
		here we assume m == us.size()-1
		*/
		constexpr int sample_num = 8;
		// u
		{
			knots_u = f_uniform_knots(*us, p, m);
			Bspline ubs(m, p + 1, knots_u);
			(*us)[0] = 0;
			(*us)[m] = 1;
			for (int i = 1; i < m; i++) {
				float delta = (knots_u[i + p + 1] - knots_u[i]) / sample_num;
				float begin = knots_u[i];
				float max_value = -1;
				float max_u = begin;
				for (int s = 0; s < sample_num; s++) {
					float value = ubs.Nik(i, p + 1, begin);
					if (value > max_value) {
						max_value = value;
						max_u = begin;
					}
					begin += delta;
				}
				(*us)[i] = max_u;
			}
		}
		// v
		{
			knots_v = f_uniform_knots(*vs, q, n);
			Bspline vbs(n, q + 1, knots_v);
			(*vs)[0] = 0;
			(*vs)[n] = 1;
			for (int i = 1; i < n; i++) {
				float delta = (knots_v[i + q + 1] - knots_v[i]) / sample_num;
				float begin = knots_v[i];
				float max_value = -1;
				float max_v = begin;
				for (int s = 0; s < sample_num; s++) {
					float value = vbs.Nik(i, q + 1, begin);
					if (value > max_value) {
						max_value = value;
						max_v = begin;
					}
					begin += delta;
				}
				(*vs)[i] = max_v;
			}
		}
	}
	// p q is degree but need order
	Bspline_Surface* bs =
		new Bspline_Surface(m, n, p + 1, q + 1, knots_u, knots_v);
	return bs;
}
