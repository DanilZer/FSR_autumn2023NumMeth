function Zernov_Task1
  clc;clear;
  Y_Data_1 = [207, 242, 271, 300, 316, 322, 320, 320, 358, 400,428, 431, 424, 406, 382, 382,  424, 420, 426, 437];
  X_Data_1 = [73, 97, 141, 194, 248, 305, 371, 420, 460, 488, 572, 628, 707, 800, 895, 997, 1039, 1052, 1098, 1123];
  figure(1); %Построение полинома Лагранжа по точкам интерполяции
  X = linspace(min(X_Data_1),max(X_Data_1));
  Y1 = lagrange(X_Data_1,Y_Data_1,X);
  plot(X_Data_1,Y_Data_1,'o',X,Y1);
  axis([0 1200 0 1000]);

  figure(2); %Построение сплайна
  [s0,s1,s2,s3] = cubic_spline(X_Data_1', Y_Data_1');
  plot(X_Data_1,Y_Data_1,"o");
  hold on;
  plot_cubic_spline(X_Data_1,s0,s1,s2,s3)
  axis([0 1200 0 1000]);

  figure(3); % построение первой производной сплайна
  hold on;
  plot_cubic_spline1d(X_Data_1,s1,s2,s3);

  figure(4); % построение второй производной сплайна
  hold on;
  plot_cubic_spline2d(X_Data_1,s2,s3);

  figure(5); % построение полного контура машины
  hold on;
  Y_FULLDATA = [207, 242, 271, 300, 316, 322, 320, 320, 358, 400,428, 431, 424, 406, 382, 382,  424, 420, 426, 437, 400, 371, 342, 318, 298, 285, 265, 252, 225, 203, 185, 172, 154, 139, 119, 101, 93, 86, 86, 101, 121, 141, 128, 124, 124, 112, 115, 115, 112, 108, 106, 101, 124, 130, 130, 119, 106, 95, 86, 79, 77, 82, 95, 110, 124, 121, 119, 99, 99, 104, 115, 132, 150, 168, 183, 188, 199,207];
  X_FULLDATA = [73, 97, 141, 194, 248, 305, 371, 420, 460, 488, 572, 628, 707, 800, 895, 997, 1039, 1052, 1098, 1123, 1107, 1089, 1083, 1081, 1081, 1096, 1092, 1085, 1045, 1012, 983, 979, 975, 968, 952, 930, 913, 882, 853, 813, 789, 771, 765, 743, 725, 705, 678, 650, 594, 517, 435, 385, 385, 373, 360, 354, 340, 327, 314, 296, 270, 236, 208, 190, 179, 161, 150, 155, 122, 102, 88, 82, 82, 73, 86, 86, 73,73];
  for (i = 1:78)
    t(i) = i;
  end
  [s0x, s1x, s2x, s3x] = cubic_spline(t', X_FULLDATA');
  [s0y, s1y, s2y, s3y] = cubic_spline(t', Y_FULLDATA');
  plot(X_FULLDATA,Y_FULLDATA,'o');
  plot_cubic_spline_xty(t, s0x, s1x, s2x, s3x, s0y, s1y, s2y, s3y)

  axis([0 1200 0 1000]);
endfunction

% функция построения интерполяционного полинома Лагранжа
function yy=lagrange(x,y,xx)
  n=length(x);
  yy=zeros(size(xx));
  for k=1:n
    t=ones(size(xx));
    for j=[1:k-1, k+1:n]
        t=t.*(xx-x(j))/(x(k)-x(j));
    end
    yy = yy + y(k)*t;
  end
endfunction



% функция для сплайна с граничными условиями
function [s0,s1,s2,s3]=cubic_spline(x,y)
% [s0,s1,s2,s3]=cubic_spline(x,y)
% Sk(x)  = sk0 + sk1*(x-x(k)) + sk2*(x-x(k))^2 + sk3*(x-x(k))^3
if any(size(x) ~= size(y)) || size(x,2) ~= 1
   error('inputs x and y must be column vectors of equal length');
end

n = length(x);

h = x(2:n) - x(1:n-1);
d = (y(2:n) - y(1:n-1))./h;

lower = h(2:end);
main  = 2*(h(1:end-1) + h(2:end));
upper = h(1:end-1);

T = spdiags([lower main upper], [-1 0 1], n-2, n-2);
T(1,1) = (3.0/2.0)*h(1)+2*h(2);%первое условие
T(end,end)=T(end,end)+h(end); %второе условие
rhs = 6*(d(2:end)-d(1:end-1));
rhs(1) = rhs(1) - 3*d(1); %первое условие

m = T\rhs;

m = [0; m; 0];
m(1) = (3.0/(h(1)))*(d(1)) - m(2)/2; %первое условие
m(n) = m(n-1); %второе условие

s0 = y; s0 = s0(1:end-1)';
s1 = d - h.*(2*m(1:end-1) + m(2:end))/6; s1 = s1(1:end)';
s2 = m/2; s2 = s2(1:end-1)';
s3 =(m(2:end)-m(1:end-1))./(6*h); s3 = s3(1:end)';
end


% функция для рисования графика
function plot_cubic_spline(x,s0,s1,s2,s3)
  n = length(x);
  inner_points = 20;
  for i=1:n-1
    xx = linspace(x(i),x(i+1),inner_points);
    xi = repmat(x(i),1,inner_points);
    yy = s0(i) + s1(i)*(xx-xi) + s2(i)*(xx-xi).^2 + s3(i)*(xx - xi).^3;
    plot(xx,yy)
  end
end


%функция для рисования графика первой производной сплайна
function plot_cubic_spline1d(x,s1,s2,s3)
  n = length(x);
  inner_points = 20;
  for i=1:n-1
    xx = linspace(x(i),x(i+1),inner_points);
    xi = repmat(x(i),1,inner_points);
    yy = s1(i) + 2*s2(i)*(xx-xi) + 3*s3(i)*(xx - xi).^2;
    plot(xx,yy)
  end
end

%функция для рисования графика второй производной сплайна
function plot_cubic_spline2d(x,s2,s3)
  n = length(x);
  inner_points = 20;
  for i=1:n-1
    xx = linspace(x(i),x(i+1),inner_points);
    xi = repmat(x(i),1,inner_points);
    yy = 2*s2(i) + 2*3*s3(i)*(xx - xi);
    plot(xx,yy)
  end
end

% функция для рисования графика сплайна заданного параметрически
function plot_cubic_spline_xty(t, s0x, s1x, s2x, s3x, s0y, s1y, s2y, s3y)
  n = length(t);
  inner_points = 20;
  for i = 1:n - 1
     tS = linspace(t(i), t(i+1), inner_points);
     t_i = repmat(t(i), 1, inner_points);
     xS = s0x(i) + s1x(i) * (tS - t_i) + s2x(i) * (tS - t_i) .^ 2 + s3x(i) * (tS - t_i) .^ 3;
     yS = s0y(i) + s1y(i) * (tS - t_i) + s2y(i) * (tS - t_i) .^ 2 + s3y(i) * (tS - t_i) .^ 3;
     plot(xS, yS)
  end
end

