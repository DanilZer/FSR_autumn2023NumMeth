function ZernovDATask2
  clc;clear;
  tspan = [0,20];
  hold on;
  figure(1);
  [B,E,vx,vy] = coef(20,5e5,1e5,1e4);
  y0 = [0;0; vx; vy];
  [t,y] = ode45(@(t,y) systemODE(B,E,t,y),tspan,y0);
  plot(y(:,1),y(:,2));
  xlabel('x');
  ylabel('y');
  title('Общий случай движения трохоида');
  figure(2);
  [B,E,vx,vy] = coef(20,5e5,0,0);
  y0 = [0;0; vx; vy];
  [t,y] = ode45(@(t,y) systemODE(B,E,t,y),tspan,y0);
  plot(y(:,1),y(:,2));
  xlabel('x');
  ylabel('y');
  title('циклоида');
  figure(3);
  [B,E,vx,vy] = coef(20,5e5,4e4,4e4);
  [B,E,vx1,vy1] = coef(20,5e5,8e4,8e4);
  [B,E,vx2,vy2] = coef(20,5e5,10e4,10e4);
  y0 = [0;0; vx; vy];
  y1 = [0;0; vx1; vy1];
  y2 = [0;0; vx2; vy2];
  [t,y] = ode45(@(t,y) systemODE(B,E,t,y),tspan,y0);
  [t,y1] = ode45(@(t,y) systemODE(B,E,t,y),tspan,y1);
  [t,y2] = ode45(@(t,y) systemODE(B,E,t,y),tspan,y2);
  plot(y(:,1),y(:,2),y1(:,1),y1(:,2),y2(:,1),y2(:,2));
  xlabel('x');
  ylabel('y');
  title('Зависимость траектории от начальной скорости');
  figure(4);
  [B,E,vx,vy] = coef(20,5e5,1e5,1e4);
  y0 = [0;0; vx; vy];
  [t,y] = ode45(@(t,y) systemODE(B,E,t,y),tspan,y0);
  Bv = [0,0,B];
  Ev = [0,E,0];
  ans = cross(Bv,Ev)/norm(B).^2;
  abs(ans(1))
  [N,ymax1] = findm(y(:,1),y(:,2),0,1);
  y(N,1)/t(N)
  [B1,E1,vx,vy] = coef(25,5e5,1e5,1e4);
  y0 = [0;0; vx; vy];
  [t,y1] = ode45(@(t,y) systemODE(B1,E1,t,y),tspan,y0);
  Bv = [0,0,B1];
  Ev = [0,E1,0];
  ans = cross(Bv,Ev)/norm(B1).^2;
  abs(ans(1))
  [N,ymax1] = findm(y1(:,1),y1(:,2),0,1);
  y1(N,1)/t(N)
  [B2,E2,vx,vy] = coef(40,5e5,1e5,1e4);
  y0 = [0;0; vx; vy];
  [t,y2] = ode45(@(t,y) systemODE(B2,E2,t,y),tspan,y0);
  Bv = [0,0,B2];
  Ev = [0,E2,0];
  ans = cross(Bv,Ev)/norm(B2).^2;
  abs(ans(1))
  [N,ymax1] = findm(y2(:,1),y2(:,2),0,1);
  y2(N,1)/t(N)
  plot(y(:,1),y(:,2),y1(:,1),y1(:,2),y2(:,1),y2(:,2));
  xlabel('x');
  ylabel('y');
  title('Зависимость траектории от скорости дрейфа');
  figure(5);
  [B,E,vx,vy] = coef(20,5e5,1e5,1e4);
  y0 = [0;0; vx; vy];
  [t,y] = ode45(@(t,y) systemODE(B,E,t,y),tspan,y0);
  [B2,E2,vx,vy] = coef(40,5e5,1e5,1e4);
  y0 = [0;0; vx; vy];
  [t2,y2] = ode45(@(t,y) systemODE(B2,E2,t,y),tspan,y0);
  plot3(t,y(:,1),y(:,2),t2,y2(:,1),y2(:,2));
  xlabel('t');
  ylabel('x');
  zlabel('y');
  title('Сравнение двух траекторий отличающихся скоростью дрейфа примерно в два раза');
endfunction

function [B,E,vx,vy] = coef(B0,E0,v1,v2)
  q = -1.602e-19;
  m = 9.31e-31;
  C = 3e8;
  l = 1e-4;
  B = (B0*q*l)/(m*C);
  E = (E0*q*l)/(m*C*C);
  vx = v1/C;
  vy = v2/C;
endfunction

function matrix = systemODE(B,E,t,y)
  A = [0,0,1,0;
      0,0,0,1;
      0,0,0,B;
      0,0,-B,0];
   matrix = A*y + [0 ; 0; 0 ; E];
endfunction

function [xmax,ymax] = findm(x,y,a,b)
  ymax = -10000;
  xmax = 0;
  for i = 1:length(x)
    if (x(i)>a && x(i)<b)
      if (y(i) > ymax)
        ymax = y(i);
        xmax = i;
      endif
    endif
  endfor
endfunction
