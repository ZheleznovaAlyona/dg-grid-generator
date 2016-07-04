#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <vector>
using namespace std;

struct Point
{
	double x;
	double y;
	bool operator==(Point point)
	{
		if(point.x == x && point.y == y)
			return true;
		else
			return false;
	}
};

struct Area
{
	int left_x;
	int right_x;
	int low_y;
	int up_y;
	bool real;
	void input(FILE *f_in);
	int ku[4]; //лево,право, низ, верх (1,2,3,-1 обозначают краевые условия 1,2, 3, никакие)
	//здесь ещё порядок базиса ввести
};

struct AreasLines
{
	vector <double> xi;
	vector <double> yj;
	void input(FILE *f_in);
};


struct Element
{
	int nodes[4];
	int edges[4]; //левое,правое, нижнее, верхнее
	int number_of_area;
	int neighbors[4]; //левый, правый, нижний, верхний
	void output(FILE *f_out);

	Element& operator=(Element element)
	{
		for(int i = 0; i < 4; i++)
			nodes[i] = element.nodes[i];
		number_of_area = element.number_of_area;
		for(int i = 0; i < 4; i++)
		neighbors[i] = element.neighbors[i];

		return *this;
	}
};

struct BoundaryCondition
{
	int elem;
	int edges[4]; //левое,правое,нижнее, верхнее: 1 - есть, 0 - нет
	int formula_number;
};

struct Partition
{
	vector <Area> areas;
	AreasLines areas_lines;
	vector <Element> elements;
	vector <Point> nodes;
	vector <BoundaryCondition> ku[3]; 

	void input(FILE *areas_f_in, FILE *areas_lines_f_in);

	void push_node(double x, double y);
	void partition_one_coordinate(vector <double> &ci,
								  vector <double> areas_lines_ci, 
								  vector <double> coefficient, 
								  vector <int> n_intervals);
	int find_area(double x, double y);
	void build_grid(FILE *intervals_f_in, vector <double> &xii, vector <double> &yjj);
	void find_and_push_neighbors(int element_number);
	void compute_elements(vector <double> xii, vector <double> yjj);
	void do_partition(FILE *intervals_f_in);
	void form_ku();
	
	void output(FILE *grid_f_out,
				FILE *elements_f_out, 
				FILE *l1_out, 
				FILE *l2_out, 
				FILE *l3_out);
};

