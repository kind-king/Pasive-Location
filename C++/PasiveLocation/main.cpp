//#include "scr/coordinates_system.h"
#include "src/Pasive_Location.h"
//#include <stdlib.h>
//#include <stdio.h>
//#include <stdbool.h>
//#include <math.h> // pi , PI

// for sample objects and output
//#include <string>
//#include <vector>
//#include <iterator> // заголовочный файл итераторов
#include <iostream>
//#include "c:\Users\IRA\documents\visual studio 2015\Projects\Sources\coordinates_system.h"
//#include <queue>
using namespace std;

int main() {
	//vector<int> myVector(10);   // объ€вл€ем вектор размером в 10 элементов и инициализируем их нул€ми
	//myVector.reserve(10);   // выдел€ем пам€ть под 10 элементов
	//myVector.push_back(25);
	//myVector.insert(myVector.end(), 3);
	//myVector.size();
	// вывод на экран элементов вектора
	//copy(myVector.begin(),   // итератор начала массива
	//	myVector.end(),     // итератор конца массива
	//	ostream_iterator<int>(cout, " ")   //итератор потока вывода
	//);
	//vector<bearing> myGoals = { { 20.0 , 1.0 , 1.0 },{ 22.0 , 1.0 , 1.0 } };
	//std::cout << "mv_1 = a: " << endl;
	//for (bearing n : myGoals) {
	//	std::cout << n.d << ' ' << n.ph << ' ' << n.z << '\n';
	//}
	//bearing bear;
	//bearing* a;
	//bearing base[4];
	bearing base[4] = { { 00.0 , 010.0 * PI / 180 ,  0.0 },
	                    { 10.0 , 200.0 * PI / 180 , -0.01 },
	                    { 10.0 , 130.0 * PI / 180 ,  0.0 },
	                    { 10.0 , 270.1 * PI / 180 ,  0.01 } };
	bearing goal;
	bearing* goal_test;
	//a[0] = bearing();
	/*for (int i = 0; i < 4; ++i) {
		std::cout << "x = " << endl;
		std::cin >> bear.d;
		std::cout << "y = " << endl;
		std::cin >> bear.ph;
		std::cout << "z = " << endl;
		std::cin >> bear.z;
		base[i] = bear;
	}*/
	Pasive_Location base1(base);
	bool exit_ = true;
	while ( exit_ ) {
		std::cout << "d = ";
		std::cin >> goal.d;
		std::cout << "ph = ";
		std::cin >> goal.ph ;
		goal.ph = goal.ph * PI / 180;
		std::cout << "z = ";
		std::cin >> goal.z;
		//bearing goal = bearing( 20.0 , 1.0 , 1.0 );
		goal_test = base1.Solve_Simulation(goal, 0);
		std::cout << "d = " << goal_test[0].d << "\tph = " << goal_test[0].ph * 180 / PI << "\tz = " << goal_test[0].z << endl;
		std::cout << "d = " << goal_test[1].d << "\tph = " << goal_test[1].ph * 180 / PI << "\tz = " << goal_test[1].z << endl;
		
	}
	
	//long double a = sqrt(-10);
	std::cout << 12 << endl;
	//std::cin >> goal.d;
	system("pause");
	return 0;
}