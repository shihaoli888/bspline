#include "fit.hpp"
#include <iostream>
int main() {
	Fitter fitter;
	fitter.knot_generation = Fitter::uniform;
	float rdata[] = {
		0.0, 0.0, 0.0,
		1.4, 0.0, 0.0,
		3.7, 0.0, 2.0,
		4.5, 0.0, 2.2,
		6.2, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.5, 1.0, 0.0,
		3.2, 1.0, 2.2,
		4.6, 1.0, 3.1,
		6.5, 1.0, 0.0,
		0.0, 2.0, 0.0,
		0.7, 2.0, 0.0,
		2.5, 2.0, 2.1,
		4.5, 2.0, 4.2,
		7.2, 2.0, 0.0,
		0.0, 3.0, 0.0,
		1.0, 3.0, 0.0,
		4.0, 3.0, 2.0,
		4.5, 3.0, 4.2,
		5.0, 3.0, 0.0 };
	int row = 4, col = 5;
	std::vector<std::vector<vec3>> data;
	data.resize(row);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			int idx = (i * col + j) * 3;
			data[i].emplace_back(rdata[idx + 0], rdata[idx + 1], rdata[idx + 2]);
		}
	}
	//int m = row - 1, n = col - 1;
	int m = row-2, n = col-2;
	std::vector<std::vector<vec3>> ps(m+1, std::vector<vec3>(n+1, vec3(0)));
	std::vector<float> us(row);
	std::vector<float> vs(col);
	Bspline_Surface* bsp = nullptr;
	//fitter.interplation(data, 2, 2, row - 1, col - 1, &ps, &us, &vs, &bsp);
	fitter.approximation(data, 2, 2, m, n, &ps, &us, &vs, &bsp);
	if (bsp) {
		auto&& bs = *bsp;
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				auto&& v = bs(ps, us[i], vs[j]);
				std::cout << v.x << " " << v.y << " " << v.z << std::endl;
			}
		}
	}
}
