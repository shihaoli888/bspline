#pragma once
#include "bspline.hpp"
struct Fitter {
	bool interplation(
		const std::vector<std::vector<vec3>>& data,
		int p,
		int q,
		int m,
		int n,
		std::vector<std::vector<vec3>>* ps,
		std::vector<float>* us,
		std::vector<float>* vs,
		Bspline_Surface** bs
	);

	bool approximation(
		const std::vector<std::vector<vec3>>& data,
		int p,
		int q,
		int m,
		int n,
		std::vector<std::vector<vec3>>* ps,
		std::vector<float>* us,
		std::vector<float>* vs,
		Bspline_Surface** bs
	);
	enum Parameter_Method
	{
		uniformly_space,
		chordlength,
		centripetal,
		universal,
	};
	enum Knot_Generation
	{
		uniform,
		average,
	};
	Parameter_Method parameter_method = uniformly_space;
	Knot_Generation knot_generation = average;
	float alpha = 0.5f;
private:
	Bspline_Surface* get_parameter(
		const std::vector<std::vector<vec3>>& data,
		int p,
		int q,
		int m,
		int n,
		std::vector<float>* us,
		std::vector<float>* vs
	)const;
};