#include "grid.h"

void Area::input(FILE *f_in)
{
	int tmp;
	fscanf(f_in, "%d", &left_x);
	fscanf(f_in, "%d", &right_x);
	fscanf(f_in, "%d", &low_y);
	fscanf(f_in, "%d", &up_y);
	fscanf(f_in, "%d", &tmp);
	if(tmp == 1) real = true;
	else real = false;
	fscanf(f_in, "%d", &ku[0]);
	fscanf(f_in, "%d", &ku[1]);
	fscanf(f_in, "%d", &ku[2]);
	fscanf(f_in, "%d", &ku[3]);
}

void AreasLines::input(FILE *f_in)
{
	int n;
	double tmp;

	fscanf(f_in, "%d", &n);
	xi.reserve(n);
	for(int i = 0; i < n; i++)
	{
		fscanf(f_in, "%lf", &tmp);
		xi.push_back(tmp);
	}

	fscanf(f_in, "%d", &n);
	yj.reserve(n);
	for(int j = 0; j < n; j++)
	{
		fscanf(f_in, "%lf", &tmp);
		yj.push_back(tmp);
	}
}

void Element::output(FILE *f_out)
{
	fprintf(f_out, "%d\n", number_of_area);

	for(int i = 0; i < 4; i++)
		fprintf(f_out, "%d ", nodes[i]);
	fprintf(f_out, "\n");
	
	for(int i = 0; i < 4; i++)
		fprintf(f_out, "%d ", neighbors[i]);
	fprintf(f_out, "\n");

	for(int i = 0; i < 4; i++)
		fprintf(f_out, "%d ", edges[i]);
	fprintf(f_out, "\n");
}

void Partition::input(FILE *areas_f_in, FILE *areas_lines_f_in)
{
	int n_areas;
	Area tmp;

	//ввод подобластей
	fscanf(areas_f_in, "%d", &n_areas);
	areas.reserve(n_areas);
	for(int i = 0; i < n_areas; i++)
	{
		tmp.input(areas_f_in);
		areas.push_back(tmp);
	}

	//ввод координатных линий
	areas_lines.input(areas_lines_f_in);
}

void Partition::push_node(double x, double y)
{
	Point p;
	p.y = y;
	p.x = x; 
	nodes.push_back(p);
}

void Partition::partition_one_coordinate(vector <double> &ci, 
										 vector <double> areas_lines_ci, 
										 vector <double> coefficient, 
										 vector <int> n_intervals)
{
	int count;
	double l; //длина интервала
	double h; //шаг
	int n_lines = areas_lines_ci.size();

	count = 0;
	for(int i = 0; i < n_lines - 1; i++)
	{
		ci.push_back(areas_lines_ci[i]);
		count++;
		
		//длина интервала
		l = abs(areas_lines_ci[i+1] - areas_lines_ci[i]);
		
		//рассчитываем первый шаг
		if(1.0 - coefficient[i] < 1E-14) 
			h = l / n_intervals[i];
		else 
			h = l * (1 - coefficient[i]) / (1 - pow(coefficient[i], n_intervals[i]));
		
		ci.push_back(ci[count - 1] + h);
		count++;
		//получаем сетку внутри интервала
		for(int j = 2; j < n_intervals[i]; j++)
		{
			h *= coefficient[i];	
			ci.push_back(ci[count - 1] + h);
			count++;
		}
	}
	ci.push_back(areas_lines_ci[n_lines - 1]);
}

int Partition::find_area(double x, double y)
{
	bool get_into_x_interval, get_into_y_interval;
	for(int i = 0; i < areas.size(); i++)
	{
		get_into_x_interval = x < areas_lines.xi[areas[i].right_x] && 
							  x > areas_lines.xi[areas[i].left_x];
		get_into_y_interval = y < areas_lines.yj[areas[i].up_y] && 
							  y > areas_lines.yj[areas[i].low_y];
		if(get_into_x_interval && get_into_y_interval) return i;
	}
}

void Partition::build_grid(FILE *intervals_f_in, vector <double> &xii, vector <double> &yjj)
{
	vector <double> coefficient_x;
	vector <double> coefficient_y;
	vector <int> n_intervals_x;
	vector <int> n_intervals_y;

	int n_lines_x = areas_lines.xi.size(), n_lines_y = areas_lines.yj.size();

	vector <double> xi; //геометрические линии разбиения по х
	vector <double> yj; //геометрические линии разбиения по у
	int nx, ny;

	int tmp1; //временная переменная для считывания
	double tmp2; //временная переменная для считывания

	int n_nodes, n_elements; //число узлов и число элементов

	n_intervals_x.reserve(n_lines_x - 1); coefficient_x.reserve(n_lines_x - 1);
	n_intervals_y.reserve(n_lines_y - 1); coefficient_y.reserve(n_lines_y - 1);
	
	//ввод количества интервалов по х и у и коэффициентов разрядки
	for(int i = 0; i < n_lines_x - 1; i++)
	{
		fscanf(intervals_f_in, "%d", &tmp1);
		n_intervals_x.push_back(tmp1);
		fscanf(intervals_f_in, "%lf", &tmp2);
		coefficient_x.push_back(tmp2);
	}

	for(int j = 0; j < n_lines_y - 1; j++)
	{
		fscanf(intervals_f_in, "%d", &tmp1);
		n_intervals_y.push_back(tmp1);
		fscanf(intervals_f_in, "%lf", &tmp2);
		coefficient_y.push_back(tmp2);
	}

	nx = 0;
	for(int i = 0; i < n_lines_x - 1; i++)
		nx += n_intervals_x[i]; //общее количество интервалов по х
	nx++;

	ny = 0;
	for(int j = 0; j < n_lines_y - 1; j++)
		ny += n_intervals_y[j];	//общее количество интервалов по у
	ny++;

	xi.reserve(nx); yj.reserve(ny);

	//построение сеток по х и у
	partition_one_coordinate(xi, areas_lines.xi, coefficient_x, n_intervals_x);
	partition_one_coordinate(yj, areas_lines.yj, coefficient_y, n_intervals_y);

	n_intervals_x.clear(); coefficient_x.clear();
	n_intervals_y.clear(); coefficient_y.clear();

	xii.reserve((nx - 1) * 3); yjj.reserve((ny - 1) * 3);

	//дублируем нужные узлы
	for(int i = 0; i < nx - 1; i++)
	{
		xii.push_back(xi[i]);
		xii.push_back(xi[i + 1]);
	}

	for(int j = 0; j < ny - 1; j++)
	{
		yjj.push_back(yj[j]);
		yjj.push_back(yj[j + 1]);
	}

	n_elements = (nx - 1) * (ny - 1);
	n_nodes = n_elements * 4;

	elements.reserve(n_elements);
	nodes.reserve(n_nodes);

	//заполняем список узлов
	for(int j = 0; j < yjj.size(); j++)
		for(int i = 0; i < xii.size(); i++)
			push_node(xii[i], yjj[j]);
}

void Partition::find_and_push_neighbors(int element_number)
{
	int count = 0, n_elements = elements.size();
	Element element;
	element = elements[element_number];
	bool tmp[4] = {false, false, false, false};

	for(int i = 0; count < 4 && i < element_number; i++)
	{
		//сосед по левому ребру
		if(nodes[element.nodes[0]] == nodes[elements[i].nodes[1]] && 
		   nodes[element.nodes[2]] == nodes[elements[i].nodes[3]])
		{
			element.neighbors[0] = i;
			tmp[0] = true;
			count++;
		}
		//сосед по правому ребру
		if(nodes[element.nodes[1]] == nodes[elements[i].nodes[0]] &&
		   nodes[element.nodes[3]] == nodes[elements[i].nodes[2]])
		{
			element.neighbors[1] = i;
			tmp[1] = true;
			count++;
		}
		//сосед по нижнему ребру
		if(nodes[element.nodes[0]] == nodes[elements[i].nodes[2]] &&
		   nodes[element.nodes[1]] == nodes[elements[i].nodes[3]])
		{
			element.neighbors[2] = i;
			tmp[2] = true;
			count++;
		}
		//сосед по верхнему ребру
		if(nodes[element.nodes[2]] == nodes[elements[i].nodes[0]] &&
		   nodes[element.nodes[3]] == nodes[elements[i].nodes[1]])
		{
			element.neighbors[3] = i;
			tmp[3] = true;
			count++;
		}
	}

	for(int i = element_number + 1; count < 4 && i < n_elements; i++)
	{
		//сосед по левому ребру
		if(nodes[element.nodes[0]] == nodes[elements[i].nodes[1]] &&
		   nodes[element.nodes[2]] == nodes[elements[i].nodes[3]])
		{
			element.neighbors[0] = i;
			tmp[0] = true;
			count++;
		}
		//сосед по правому ребру
		if(nodes[element.nodes[1]] == nodes[elements[i].nodes[0]] &&
		   nodes[element.nodes[3]] == nodes[elements[i].nodes[2]])
		{
			element.neighbors[1] = i;
			tmp[1] = true;
			count++;
		}
		//сосед по нижнему ребру
		if(nodes[element.nodes[0]] == nodes[elements[i].nodes[2]] &&
		   nodes[element.nodes[1]] == nodes[elements[i].nodes[3]])
		{
			element.neighbors[2] = i;
			tmp[2] = true;
			count++;
		}
		//сосед по верхнему ребру
		if(nodes[element.nodes[2]] == nodes[elements[i].nodes[0]] &&
		   nodes[element.nodes[3]] == nodes[elements[i].nodes[1]])
		{
			element.neighbors[3] = i;
			tmp[3] = true;
			count++;
		}
	}

	for(int i = 0; i < 4; i++)
		if(tmp[i] == false) element.neighbors[i] = -1; //отсутствие соседа на соответствующей стороне

	elements[element_number] = element;
}

void Partition::compute_elements(vector <double> xii, vector <double> yjj)
{
	int n_intervals_x = xii.size() / 2, n_intervals_y = yjj.size() / 2;
	int ii, jj;
	Element tmp_el;
	double x, y;

	//вычисляем глобальные номера узлов каждого кэ и заполняем список кэ
	for(int j = 0; j < n_intervals_y; j++)
	{
		for(int i = 0; i < n_intervals_x; i++)
		{
			for(int iy = 0; iy < 2; iy++)
			{
				jj = 2 * j + iy;
				for(int ix = 0; ix < 2; ix++)
				{					
					ii = 2 * i +ix;
					tmp_el.nodes[iy * 2 + ix] = 2 * jj * n_intervals_x + ii;
				}
			}
			int low = j * n_intervals_x * 4 + i;
			int left = j * n_intervals_x * 4 + n_intervals_x + i * 2;
			int right = left + 1;
			int up = j * n_intervals_x * 4 + n_intervals_x + n_intervals_x * 2 + i;
			tmp_el.edges[0] = left;
			tmp_el.edges[1] = right;
			tmp_el.edges[2] = low;
			tmp_el.edges[3] = up;
			elements.push_back(tmp_el);
		}
	}

	//находим номер подобласти для каждого кэ
	for(int i = 0; i < elements.size(); i++)
	{
		//центральный узел берём
		x = (nodes[elements[i].nodes[0]].x + nodes[elements[i].nodes[1]].x) / 2;
		y = (nodes[elements[i].nodes[0]].y + nodes[elements[i].nodes[2]].y) / 2;
		elements[i].number_of_area = find_area(x, y);		
	}

	//находим соседей
	for(int i = 0; i < elements.size(); i++)
		find_and_push_neighbors(i);
}

void Partition::do_partition(FILE *intervals_f_in)
{
	vector <double> xii; //линии разбиения с учётом дублирования узлов по х
	vector <double> yjj; //линии разбиения с учётом дублирования узлов по у
	build_grid(intervals_f_in, xii, yjj);
	compute_elements(xii, yjj);
}

void Partition::form_ku()
{
	int size = elements.size();
	BoundaryCondition tmp;

	//для каждого элемента находим,какие ку на его границах
	for(int i = 0; i < size; i++)
	{
		bool left[3] = {false, false, false}, right[3] = {false, false, false}, 
			 low[3] = {false, false, false}, up[3] = {false, false, false};
		if(elements[i].neighbors[0] == -1)
		{
			if(areas[elements[i].number_of_area].ku[0] == 1) left[0] = true;
			else
				if(areas[elements[i].number_of_area].ku[0] == 2) left[1] = true;
				else
					if(areas[elements[i].number_of_area].ku[0] == 3) left[2] = true;
		}
		if(elements[i].neighbors[1] == -1)
		{
			if(areas[elements[i].number_of_area].ku[1] == 1) right[0] = true;
			else
				if(areas[elements[i].number_of_area].ku[1] == 2) right[1] = true;
				else
					if(areas[elements[i].number_of_area].ku[1] == 3) right[2] = true;
		}
		if(elements[i].neighbors[2] == -1)
		{
			if(areas[elements[i].number_of_area].ku[2] == 1) low[0] = true;
			else
				if(areas[elements[i].number_of_area].ku[2] == 2) low[1] = true;
				else
					if(areas[elements[i].number_of_area].ku[2] == 3) low[2] = true;
		}
		if(elements[i].neighbors[3] == -1)
		{
			if(areas[elements[i].number_of_area].ku[3] == 1) up[0] = true;
			else
				if(areas[elements[i].number_of_area].ku[3] == 2) up[1] = true;
				else
					if(areas[elements[i].number_of_area].ku[3] == 3) up[2] = true;
		}

		for(int j = 0; j < 3; j++)
			if(left[j] || right[j] || low[j] || up[j])
			{
				tmp.elem = i;
				tmp.formula_number = elements[i].number_of_area;
				if(left[0]) tmp.edges[0] = 1;
				else tmp.edges[0] = 0;
				if(right[0]) tmp.edges[1] = 1;
				else tmp.edges[1] = 0;
				if(low[0]) tmp.edges[2] = 1;
				else tmp.edges[2] = 0;
				if(up[0]) tmp.edges[3] = 1;
				else tmp.edges[3] = 0;
				ku[j].push_back(tmp);
			}
	}
}

void Partition::output(FILE *grid_f_out, 
					   FILE *elements_f_out, 
					   FILE *l1_out, 
					   FILE *l2_out, 
					   FILE *l3_out)
{
	int n_nodes = nodes.size();
	int n_l1 = ku[0].size(), n_l2 = ku[1].size(), n_l3 = ku[2].size();

	//выводим список узлов
	fprintf(grid_f_out, "%d\n", n_nodes);
	for(int i = 0; i < n_nodes; i++)
		fprintf(grid_f_out, "%.20lf %.20lf\n", nodes[i].x, nodes[i].y);

	//выводим список элементов
	for(int i = 0; i < elements.size(); i++)
		elements[i].output(elements_f_out);

	fprintf(l1_out, "%d\n", n_l1);
	fprintf(l2_out, "%d\n", n_l2);
	fprintf(l3_out, "%d\n", n_l3);
	for(int i = 0; i < n_l1; i++)
	{
		fprintf(l1_out, "%d ", ku[0][i].elem);
		fprintf(l1_out, "%d\n", ku[0][i].formula_number);
		fprintf(l1_out, "%d %d %d %d\n", 
				ku[0][i].edges[0], ku[0][i].edges[1], ku[0][i].edges[2], ku[0][i].edges[3]);
	}

	for(int i = 0; i < n_l2; i++)
	{
		fprintf(l2_out, "%d ", ku[1][i].elem);
		fprintf(l2_out, "%d\n", ku[1][i].formula_number);
		fprintf(l2_out, "%d %d %d %d\n", 
				ku[1][i].edges[0], ku[1][i].edges[1], ku[1][i].edges[2], ku[1][i].edges[3]);
	}

	for(int i = 0; i < n_l3; i++)
	{
		fprintf(l3_out, "%d ", ku[2][i].elem);
		fprintf(l3_out, "%d\n", ku[2][i].formula_number);
		fprintf(l3_out, "%d %d %d %d\n", 
				ku[2][i].edges[0], ku[2][i].edges[1], ku[2][i].edges[2], ku[2][i].edges[3]);
	}

}

void main()
{
	Partition partition;
	FILE *areas_f_in, *areas_lines_f_in, *intervals_f_in;
	FILE *grid_f_out, *elements_f_out;
	FILE *l1_out, *l2_out, *l3_out;

	areas_f_in = fopen("areas.txt", "r");
	areas_lines_f_in = fopen("lines.txt", "r");
	intervals_f_in = fopen("ints.txt", "r");
	grid_f_out = fopen("grid.txt", "w"); 
	elements_f_out = fopen("elements.txt", "w");
	l1_out = fopen("l1.txt", "w");
	l2_out = fopen("l2.txt", "w");
	l3_out = fopen("l3.txt", "w");
	
	partition.input(areas_f_in, areas_lines_f_in);
	partition.do_partition(intervals_f_in);
	partition.form_ku();
	partition.output(grid_f_out, elements_f_out, l1_out, l2_out, l3_out);

	fclose(areas_f_in);
	fclose(areas_lines_f_in);
	fclose(intervals_f_in);
	fclose(grid_f_out);
	fclose(elements_f_out);
	fclose(l1_out);
	fclose(l2_out);
	fclose(l3_out);

	_getch();
}