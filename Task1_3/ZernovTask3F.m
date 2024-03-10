function ZernovTask3F
  clc; clear;
  A = 1.5;
  w0 = 2/3;
  Q = 1.3;
  A = 1.5;
  x00 = [0;0];
  x01 = [0;-3];
  T1 = 200;
  T2 = 1000;
  N = 10000;
  h = (T2 - T1) / N;
  tspan = T1 : h : T2 - h;
  [t1s,x1s] = ode45(@(t,x) system(t,x,Q,w0,A),[0,T1],x00);
  x_1s = x1s(end, :);
  [t1,x1] = ode45(@(t,x) system(t,x,Q,w0,A), tspan, x_1s);
  [t2s,x2s] = ode45(@(t,x) system(t,x,Q,w0,A),[0,T1],x01);
  x_2s = x2s(end, :);
  [t2,x2] = ode45(@(t,x) system(t,x,Q,w0,A), tspan, x_2s);
  Y1 = fft(x1(:, 1));
  P1 = 1 / N * sqrt((Y1 .* conj(Y1)));
  Y2 = fft(x2(:, 1));
  P2 = 1 / N * sqrt((Y2 .* conj(Y2)));
  w = 0 : 2 * pi / (T2 - T1) : 2 * pi / (T2 - T1) * (N - 1);
  figure(1);
  plot(x1(:,1),x1(:,2));
  xlabel("x1");
  ylabel("x2");
  title("Фазовая кривая для x(0)=(0,0)");
  figure(2);
  plot(x2(:,1),x2(:,2));
  xlabel("x1");
  ylabel("x2");
  title("Фазовая кривая для x(0)=(0,-3)");
  figure(3);
  [w1m,p1m] = findm(w,P1,0,10)
  plot(w(1:N/20),P1(1:N/20),w1m,p1m,"o");
  T1 = 2*pi/w1m
  xlabel("w");
  ylabel("p");
  title("Спектограмма для x1 для x(0)=(0,0)");
  figure(4);
  [w2m,p2m] = findm(w,P2,0,10)
  plot(w(1:N/20),P2(1:N/20),w2m,p2m,"o");
  T2 = 2*pi/w2m
  xlabel("w");
  ylabel("p");
  title("Спектограмма для x1 для x(0)=(0,-3)");
  figure(5);
  [x1m,x1y] = findm(t1,x1(:,1),905,910);
  [x2m,x2y] = findm(t1,x1(:,1),915,920);
  T1 = abs(x2m-x1m)
  plot(t1(8800:9000),x1(8800:9000,1),x1m,x1y,"o",x2m,x2y,"o");
  xlabel("t");
  ylabel("x1");
  [x1m,x1y] = findm(t1,x1(:,1),200,210);
  [x2m,x2y] = findm(t1,x1(:,1),211,220);
  T1 = abs(x2m-x1m)
  figure(6);
  [x1m,x1y] = findm(t2,x2(:,1),905,910);
  [x2m,x2y] = findm(t2,x2(:,1),915,920);
  plot(t2(8800:9000),x2(8800:9000),x1m,x1y,"o",x2m,x2y,"o");
  xlabel("t");
  ylabel("x1");
  T2 = abs(x2m-x1m)
endfunction

function xans = system(t,x,Q,w,A)
  ans1 = x(2);
  ans2 = (-1/Q)*x(2) - sin(x(1)) + A*cos(w*t);
  xans = [ans1 ; ans2];
endfunction

function [xmax,ymax] = findm(x,y,a,b)
  ymax = -10000;
  xmax = 0;
  for i = 1:length(x)
    if (x(i)>a && x(i)<b)
      if (y(i) > ymax)
        ymax = y(i);
        xmax = x(i);
      endif
    endif
  endfor
endfunction
