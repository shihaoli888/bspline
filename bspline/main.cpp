#include "fit.hpp"
#include <glad/glad.h>
#include <GL/GL.h>
#include <GLFW/glfw3.h>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>


namespace {
	constexpr int WIDTH = 800, HEIGHT = 800;
	//窗口显示相关
	GLuint vao[3];
	GLuint vbo[3];
	GLuint ebo[3];
	//GLuint vbo, ebo;
	GLuint program;
	int index_size[3];
	const char* vs = R"(#version 330
	layout (location = 0) in vec3 vPosition;
    layout (location = 1) in vec3 vNormal;
	out vec3 world_normal;
    uniform mat4 V;
    uniform mat4 P;
	void main()
	{
        world_normal = vNormal;
		gl_Position = P*V*vec4(vPosition, 1.0);
        //gl_Position = vec4(vPosition,1.0);
	}
	)";
	const char* fs = R"(#version 330
	in vec3 world_normal;
	out vec4 fragment_color;
	uniform int shading_type;
	void main(){
        if(shading_type==0){
			fragment_color = vec4(1,1,1,1);
			return;
		}
        if(shading_type==1){
			fragment_color=vec4(0,1,0,1);
			return;
		}
		vec3 n = normalize(world_normal);
		vec3 color = (n+vec3(1,1,1))/2;
		fragment_color = vec4(color,1);
	}
	)";
	//算法和模型相关
	vec3 bound_min(INFINITY);
	vec3 bound_max(-INFINITY);
	float bound_radius;
	vec3 bound_center;
	Fitter fitter;
	Bspline_Surface* bs = nullptr;
	std::vector<std::vector<vec3>>* data;
	std::vector<float>* dus, * dvs;
	std::vector<std::vector<vec3>>* ps;
	//控制相关
	float theta = 1.57;
	float phi = 0;
	bool change = false;
	bool is_approx = false;
	float t_step = 0.02f;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
	else if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
		theta -= 0.157;
		theta = std::max(0.1f, theta);
	}
	else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
		theta += 0.157;
		theta = std::min(3.13f, theta);
	}
	else if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
		phi += 0.157;
		phi = std::min(6.28f, phi);
	}
	else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
		phi -= 0.157;
		phi = std::max(0.01f, phi);
	}
	else if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
		if (fitter.parameter_method == Fitter::uniformly_space) {
			std::cout << "parameter method: uniform -> chord" << std::endl;
			fitter.parameter_method = Fitter::chordlength;
			change = true;
		}
		else if (fitter.parameter_method == Fitter::chordlength) {
			std::cout << "parameter method: chord -> centripetal 0.5" << std::endl;
			fitter.parameter_method = Fitter::centripetal;
			change = true;
		}
		else if (fitter.parameter_method == Fitter::centripetal) {
			if (!is_approx) {
				std::cout << "parameter method: centripetal 0.5 -> universal" << std::endl;
				fitter.parameter_method = Fitter::universal;
			}
			else {
				std::cout << "parameter method: centripetal 0.5 -> uniform" << std::endl;
				fitter.parameter_method = Fitter::uniformly_space;
			}
			change = true;
		}
		else if (fitter.parameter_method == Fitter::universal) {
			std::cout << "parameter method: universal -> uniform" << std::endl;
			fitter.parameter_method = Fitter::uniformly_space;
			change = true;
		}
		else;
	}
	else if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) {
		if (fitter.knot_generation == Fitter::uniform) {
			std::cout << "knots generator method: uniform -> average" << std::endl;
			fitter.knot_generation = Fitter::average;
			change = true;
		}
		else if (fitter.knot_generation == Fitter::average) {
			std::cout << "knots generator method: average -> uniform" << std::endl;
			fitter.knot_generation = Fitter::uniform;
			change = true;
		}
		else;
	}
}
void init_opengl_data() {
	//create vao vbo ebo
	glGenVertexArrays(3, vao);
	glGenBuffers(3, vbo);
	glGenBuffers(3, ebo); 
}

void refrash_opengl_data() {
	int row = ps->size();
	int col = (*ps)[0].size();
	auto&& rps = *ps;
	auto&& rbs = *bs;
	auto&& rdata = *data;
	int data_row = rdata.size();
	int data_col = rdata[0].size();
	//point
	glBindVertexArray(vao[0]);
	std::vector<float> point_data(6 * data_row * data_col);
	for (int i = 0; i < data_row; i++)
		for (int j = 0; j < data_col; j++)
			for (int k = 0; k < 3; k++)
				point_data[(i * data_col + j) * 6 + k] = rdata[i][j][k];
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, point_data.size() * sizeof(float), point_data.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	std::vector<int> point_index(data_row * data_col);
	for (int i = 0; i < data_row * data_col; i++)point_index[i] = i;
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, point_index.size() * sizeof(int), point_index.data(), GL_STATIC_DRAW);
	index_size[0] = point_index.size();
	//line
	glBindVertexArray(vao[1]);
	std::vector<float> line_data(6 * row * col);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			for (int k = 0; k < 3; k++)
				line_data[(i * col + j) * 6 + k] = rps[i][j][k];
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, line_data.size() * sizeof(float), line_data.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[1]);
	std::vector<int> line_index;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col-1; j++) {
			int idx = i * col + j;
			line_index.push_back(idx);
			line_index.push_back(idx + 1);
		}
	}
	for (int j = 0; j < col; j++) {
		for (int i = 0; i < row - 1; i++) {
			int idx = i * col + j;
			line_index.push_back(idx);
			line_index.push_back(idx + col);
		}
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, line_index.size() * sizeof(int), line_index.data(), GL_STATIC_DRAW);
	//glBindVertexArray(0);
	index_size[1] = line_index.size();
	//face
	bound_min = vec3(INFINITY);
	bound_max = vec3(-INFINITY);
	glBindVertexArray(vao[2]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	std::vector<vec3> face_data;
	for (float ut = 0; ut <= 1.f; ut += t_step) {
		for (float vt = 0; vt <= 1.f; vt += t_step) {
			auto&& p = rbs(rps, ut, vt);
			bound_min.x = std::min(bound_min.x, p.x);
			bound_min.y = std::min(bound_min.y, p.y);
			bound_min.z = std::min(bound_min.z, p.z);
			bound_max.x = std::max(bound_max.x, p.x);
			bound_max.y = std::max(bound_max.y, p.y);
			bound_max.z = std::max(bound_max.z, p.z);
			face_data.push_back(p);
			face_data.push_back(vec3(0));
		}
	}
	bound_center = (bound_max + bound_min) / 2.f;
	bound_radius = glm::length(bound_max - bound_min) / 2.f;
	int n = face_data.size() / 2;
	n = sqrt(n);
	assert(n * n * 2 == face_data.size());
	for (int r = 0; r < n-1; r ++) {
		for (int c = 0; c < n-1; c++) {
			int idx = (r * n + c) * 2;
			auto&& p0 = face_data[idx];
			auto&& p1 = face_data[idx + n * 2];
			auto&& p2 = face_data[idx + 2];
			auto&& normal = glm::normalize(glm::cross(p1 - p0, p2 - p0));
			face_data[idx + 1] = normal;
			face_data[idx + 3] = normal;
			face_data[idx + n * 2 + 1] = normal;
			face_data[idx + n * 2 + 3] = normal;
		}
	}
	glBufferData(GL_ARRAY_BUFFER, face_data.size() * sizeof(vec3), face_data.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo[2]);
	std::vector<int> face_index;
	for (int r = 0; r < n - 1; r++) {
		for (int c = 0; c < n - 1; c++) {
			int idx = r * n + c;
			face_index.push_back(idx);
			face_index.push_back(idx + n);
			face_index.push_back(idx + 1);
			face_index.push_back(idx + n);
			face_index.push_back(idx + n + 1);
			face_index.push_back(idx + 1);
		}
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, face_index.size() * sizeof(int), face_index.data(), GL_STATIC_DRAW);
	index_size[2] = face_index.size();
	glBindVertexArray(0);
}

void init_shader() {
	GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(v_shader, 1, &vs, NULL);
	glCompileShader(v_shader);
	int success = 0;
	char infoLog[512];
	glGetShaderiv(v_shader, GL_COMPILE_STATUS, &success);

	if (!success) {
		glGetShaderInfoLog(v_shader, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
	}

	GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(f_shader, 1, &fs, NULL);
	glCompileShader(f_shader);
	glGetShaderiv(f_shader, GL_COMPILE_STATUS, &success);

	if (!success) {
		glGetShaderInfoLog(f_shader, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
	}

	program = glCreateProgram();
	glAttachShader(program, v_shader);
	glAttachShader(program, f_shader);
	glLinkProgram(program);
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(program, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
	}
}

void render() {
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	float cos_theta = std::cos(theta);
	float cos_phi = std::cos(phi);
	float sin_theta = std::sin(theta);
	float sin_phi = std::sin(phi);
	vec3 dir{ sin_theta * sin_phi,cos_theta,sin_theta * cos_phi };
	vec3 eye = bound_center + 2.f * bound_radius * dir;
	float tnear = 0.5f * bound_radius, tfar = 3.5f * bound_radius;
	auto view = glm::lookAt(eye, bound_center, { 0.f, 1.f, 0.f });
	auto perspective = glm::perspective(glm::radians(60.f), WIDTH * 1.f / HEIGHT, tnear, tfar);
	glUseProgram(program);
	glUniformMatrix4fv(glGetUniformLocation(program, "P"), 1, GL_FALSE, &(perspective[0][0]));
	glUniformMatrix4fv(glGetUniformLocation(program, "V"), 1, GL_FALSE, &(view[0][0]));
	//draw face
	glUniform1i(glGetUniformLocation(program, "shading_type"), 2);
	glBindVertexArray(vao[2]);
	glDrawElements(GL_TRIANGLES, index_size[2], GL_UNSIGNED_INT, 0);
	//draw line
	glUniform1i(glGetUniformLocation(program, "shading_type"), 1);
	glBindVertexArray(vao[1]);
	glDrawElements(GL_LINES, index_size[1], GL_UNSIGNED_INT, 0);
	//draw point
	glUniform1i(glGetUniformLocation(program, "shading_type"), 0);
	glBindVertexArray(vao[0]);
	glDrawElements(GL_POINTS, index_size[0], GL_UNSIGNED_INT, 0);
}

bool read_file(const char* file_name,std::vector<std::vector<vec3>>* data, int* row, int* col) {
	std::ifstream in(file_name);
	if (!in.is_open()) return false;
	in >> (*row) >> (*col);
	data->resize(*row);
	for (int i = 0; i < *row; i++) {
		(*data)[i].resize(*col);
	}
	for (int i = 0; i < *row; i++) {
		for (int j = 0; j < *col; j++) {
			vec3 p;
			in >> p.x >> p.y >> p.z;
			(*data)[i][j] = p;
		}
	}
	return true;
}
int main(int argc, char** argv) {
	int m, n;
	int p, q;
	if (strcmp(argv[1], "interplation") == 0) {
		is_approx = false;
	}
	else if (strcmp(argv[1], "approximate") == 0) {
		is_approx = true;
	}
	else {
		std::cerr << "unknow type" << std::endl;
		return -1;
	}
	int num_row, num_col;
	data = new std::vector<std::vector<vec3>>();
	if (!read_file(argv[2], data, &num_row, &num_col)) return -1;
	p = std::stoi(argv[3]);
	q = std::stoi(argv[4]);
	if (!is_approx) {
		m = num_row - 1;
		n = num_col - 1;
	}
	else {
		m = std::stoi(argv[5]);
		n = std::stoi(argv[6]);
	}
	//init
	ps = new std::vector<std::vector<vec3>>();
	ps->resize(m + 1);
	for (int i = 0; i < m + 1; i++) (*ps)[i].resize(n + 1);
	dus = new std::vector<float>(num_row);
	dvs = new std::vector<float>(num_col);
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Bspline Surface fit", nullptr, nullptr);
	if (window == nullptr)
	{
		std::cerr << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}
	init_shader();
	init_opengl_data();
	glViewport(0, 0, WIDTH, HEIGHT);
	glEnable(GL_DEPTH_TEST);
	glPointSize(4);
	change = true;
	while (!glfwWindowShouldClose(window))
	{
		glfwPollEvents();
		if (change) {
			delete bs;
			if (!is_approx) {
				if (fitter.interplation(*data, p, q, m, n, ps, dus, dvs, &bs)) {
					refrash_opengl_data();
				}
			}
			else {
				if (fitter.approximation(*data, p, q, m, n, ps, dus, dvs, &bs)) {
					refrash_opengl_data();
				}
			}
			change = false;
		}
		render();
		glfwSwapBuffers(window);
	}
	glfwTerminate();
	delete bs;
	delete ps;
	delete dus;
	delete dvs;
	delete data;
	return 0;
}
