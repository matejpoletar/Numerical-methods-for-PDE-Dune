Point(1) = {0,5,0.1};
Point(2) = {15,5,0.1};
Point(3) = {15,-5,0.1};
Point(4) = {0,-5,0.1};
Point(5) = {0,0,0.1}; //srediste kruznice
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Circle(4) = {1,5,4};
Line Loop (1) = {1,2,3,4};

Point(6) = {0,0.5,0.1};
Circle(5) = {6,5,6};

Line Loop(6) = {4, -3, -2, -1};
Line Loop(7) = {5};
Plane Surface(8) = {6, 7};


