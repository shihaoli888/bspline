#pragma once
#include <glm/vec3.hpp>
#include <vector>
#include <assert.h>

using glm::vec3;
struct Bspline {
	Bspline(int n_, int k_, const std::vector<float>& knots_) :
		n(n_), k(k_), knots(knots_) {
		assert(n + k + 1 <= knots.size());
	}
	vec3 operator()(const std::vector<vec3>& ps, float t) const {
		//DeBoor-Cox calculation
		if (t < knots[k - 1] || t >= knots[n + 1]) return vec3(0);
		auto&& nik = Nik(k, t);
		vec3 res(0,0,0);
		for (int i = 0; i < nik.size(); i++) res += nik[i] * ps[i];
		return res;
		/*std::vector<vec3> tps(ps);
		for (int j = 1; j <= k - 1; j++) {
			for (int i = n; i >= j; i--) {
				float denominator = knots[i + k - j] - knots[i];
				if (denominator == 0) tps[i] = vec3(0);
				else {
					float tij = (t - knots[i]) / denominator;
					tps[i] = (1 - tij) * tps[i - 1] + tij * tps[i];
				}
			}
		}
		int i = k-1;
		for (; i < n + 1; i++) {
			if (t < knots[i + 1]) break;
		}
		return tps[i];*/
	}
	float Nik(int i, int k, float t) const {
		if (t < knots[i] || t >= knots[i + k]) return 0;
		std::vector<float> dp(n + 1, 0);
		int j = 0;
		for (j = 0; j <= n; j++) {
			if (t >= knots[j] && t < knots[j + 1]) {
				dp[j] = 1;
				break;
			}
		}
		for (int kk = 2; kk <= k; kk++) {
			for (int ii = 0; ii <= n; ii++) {
				float a = dp[ii];
				float b = ii == n ? 0 : dp[ii + 1];
				float denominator1 = knots[ii + k - 1] - knots[ii];
				float factor1 = t - knots[ii];
				if (denominator1 == 0) factor1 = 0;
				else factor1 /= denominator1;
				float denominator2 = knots[ii + k] - knots[ii + 1];
				float factor2 = knots[ii + k] - t;
				if (denominator2 == 0) factor2 = 0;
				else factor2 /= denominator2;
				dp[ii] = factor1 * a + factor2 * b;
			}
		}
		return dp[i];
	}
	std::vector<float> Nik(int k, float t) const {
		std::vector<float> dp(n + 1, 0);
		if (t < knots[k - 1] || t >= knots[n + 1]) return dp;
		int j = 0;
		for (j = 0; j <= n; j++) {
			if (t >= knots[j] && t < knots[j + 1]) {
				dp[j] = 1;
				break;
			}
		}
		for (int kk = 2; kk <= k; kk++) {
			for (int ii = 0; ii <= n; ii++) {
				float a = dp[ii];
				float b = ii == n ? 0 : dp[ii + 1];
				float denominator1 = knots[ii + k - 1] - knots[ii];
				float factor1 = t - knots[ii];
				if (denominator1 == 0) factor1 = 0;
				else factor1 /= denominator1;
				float denominator2 = knots[ii + k] - knots[ii + 1];
				float factor2 = knots[ii + k] - t;
				if (denominator2 == 0) factor2 = 0;
				else factor2 /= denominator2;
				dp[ii] = factor1 * a + factor2 * b;
			}
		}
		return dp;
	}
	int n, k;
	const std::vector<float> knots;
};

struct Bspline_Surface {
	Bspline_Surface(int m_, int n_, int p_, int q_,
		const std::vector<float>& knots_u_, const std::vector<float>& knots_v_)
		:ubs(m_,p_,knots_u_),vbs(n_,q_,knots_v_){
	}
	vec3 operator()(const std::vector<std::vector<vec3>>& ps, float u, float v) const {
		std::vector<vec3> tmp;
		for (auto&& p : ps) {
			tmp.emplace_back(vbs(p, v));
		}
		return ubs(tmp, u);
	}
	const Bspline ubs, vbs;
};