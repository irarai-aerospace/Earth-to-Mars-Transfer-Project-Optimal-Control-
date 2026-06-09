A = [-5, -11, 5, -6, 1 ,0, -6, 6, -6, -3;
      8 , 2 ,-5,  0 , 1 ,-2, 8, -3, 7 ,4;
      1, -6, 3, 4, -3, -4, -1 ,11, -3, -8;
      0, -4, 9, 0 ,-7, -3,  -10, -5 ,3, 0;
      8, 2, -5, -5, 3, -5, -3, -5, -10, 3];

B = [-5, -11, 5,-6, 1 ;
      8 , 2 ,-5, 0, 1 ;
      1, -6, 3, 4, -3;
      0, -4, 9, 0 ,-7;
      8, 2, -5, -5, 3];

y = [1;2;3;4;5];
%% a.1 convex

cvx_begin
    variable x(10);
    minimize(norm(x,2));
    subject to
       norm(A*x + B*y, 2) <= 0.1;
cvx_end

%% a.2  not convex

%% a.3 convex

cvx_begin
    variable x(10);
    minimize(sum_square(x));
    subject to
       norm(A*x + B*y, 2) <= 0.1;
cvx_end

%% a.4 not convex
%% a.5 not convex
%% a.6 
one = ones(10,1);
cvx_begin
    variable x(5);
    minimize(norm(A - (x + y)*one', 2));    
cvx_end

%% a.7

cvx_begin
    variable X(5,5) symmetric
    minimize(norm(X -B,'fro'));
    subject to
       X >= 0;
cvx_end

%% a.8 

cvx_begin
    variable X(5,5) symmetric
    minimize( norm(X - B,'fro'))
    subject to
        X <= 0;   
cvx_end
%% a.9 

cvx_begin
    variable X(5,5) symmetric
    minimize(norm(X - B, 2));
    subject to
       0 <= X <= eye(5);
cvx_end
