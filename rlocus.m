s = tf('s');
G = 1 / (s*(s+1)*(s+2));

% 绘制根轨迹
figure;
rlocus(G);
grid on;
title('根轨迹');

% 计算渐近线
poles = pole(G);   %开环极点
zeros = zero(G);   %开环零点

[poles_conj,zeros_conj] = is_conjugate(poles,zeros);  %提取共轭零极点
n = length(poles);
m = length(zeros);
sigma = (sum(poles) - sum(zeros)) / (n - m);
angles = (2 * (0:(n - m - 1)) + 1) * 180 / (n - m);

% 绘制渐近线
hold on;
for angle_i = angles
    x = sigma + [0, cosd(angle_i) * 10];
    y = [0, sind(angle_i) * 10];
    plot(x, y, '--r','color','#7E2F8E');
end

% 计算分离点（符号计算）
syms s K;
[num, den] = tfdata(G, 'v');
char = poly2sym(den, s) + K * poly2sym(num, s);
char_eq = poly2sym(den, s) + K * poly2sym(num, s) == 0;
dK_ds = diff(solve(char_eq, K), s);
separation_points = double(solve(dK_ds == 0, s));

% 在图中标注分离点
plot(real(separation_points), imag(separation_points), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
text(real(separation_points), imag(separation_points), ' 分离/汇合点', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
for j = 1:length(separation_points)
    text(real(separation_points(j)), imag(separation_points(j)), sprintf("%.2f+%.2fj",real(separation_points(j)), imag(separation_points(j))), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end
% 计算极点处的出射角
if ~isempty(poles_conj) 
    for k = 1:length(poles_conj)
        selected_pole = poles_conj(k);
        angle_sum_poles = 0;
        angle_sum_zeros = 0;
        % 计算到其他极点的角度之和
        for j = 1:length(poles)
            if poles(j) ~= selected_pole
                angle_temp = angle(selected_pole - poles(j)) * (180/pi);
                angle_sum_poles = angle_sum_poles + angle_temp;
            end
        end
        % 计算到所有零点的角度之和
        if ~isempty(zeros)
            for j = 1:length(zeros)
                angle_temp = angle(selected_pole - zeros(j)) * (180/pi);
                angle_sum_zeros = angle_sum_zeros + angle_temp;
            end
        else
            angle_sum_zeros = 0;
        end
        
        
        % 计算出射角
        departure_angle = 180 + angle_sum_zeros-angle_sum_poles ;
        departure_angle = mod(departure_angle + 180, 360) - 180;
        % 在图中标注出射角
        quiver(real(selected_pole), imag(selected_pole), cosd(departure_angle), sind(departure_angle), 'AutoScale', 'off', 'MaxHeadSize', 2, 'Color', 'black');
        text(real(selected_pole), imag(selected_pole), sprintf(' %.2f°', departure_angle), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end
end

if ~isempty(zeros_conj) 
    for k = 1:length(zeros_conj)
    selected_zero = zeros_conj(k);
    angle_sum_poles = 0;
    angle_sum_zeros = 0;
    % 计算到所有极点的角度之和
    for i = 1:length(poles)
        angle_temp = angle(selected_zero - poles(i)) * (180/pi);
        angle_sum_poles = angle_sum_poles + angle_temp;
    end
    % 计算到其他零点的角度之和
    for i = 1:length(zeros)
        if zeros(i) ~= selected_zero
            angle_temp = angle(selected_zero - zeros(i)) * (180/pi);
            angle_sum_zeros = angle_sum_zeros + angle_temp;
        end
    end
    % 计算入射角
    arrival_angle = 180 -( angle_sum_zeros-angle_sum_poles );
    arrival_angle = mod(arrival_angle + 180, 360) - 180;
    % 在图中标注入射角
    quiver(real(selected_zero), imag(selected_zero), cosd(arrival_angle), sind(arrival_angle), 'AutoScale', 'off', 'MaxHeadSize', 2, 'Color', 'm');
    text(real(selected_zero), imag(selected_zero), sprintf(' %.2f°', arrival_angle), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end
end

% 构建闭环特征方程
Ds = flip(coeffs(expand(char),s));
[Gain, w] = findPureImagRoots(Ds);
if ~isempty(w)
    plot(0, w, 'r*', 'MarkerSize', 10, 'DisplayName', 'jw');
    for j = 1:length(w)
        text(0, w(j), sprintf('%.2fj', w(j)), 'Color', 'red', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end
end
if ~isempty(zeros)
    x_min = min(min(real(poles)),min(real(zeros)));
    x_max = max(max(real(poles)),max(real(zeros)));
    y_min = min(min(imag(poles)),min(imag(zeros)));
    y_max = max(max(imag(poles)),max(imag(zeros)));

else
    x_min = min(real(poles));
    x_max = max(real(poles));
    y_min = min(imag(poles));
    y_max = max(imag(poles));
end
x_range = x_max - x_min;
y_range = y_max - y_min;
x_plot = [x_min - 0.25*x_range,x_max+0.25*x_range];
y_plot = [y_min - 0.25*y_range,y_max+0.25*y_range];
if(x_range == 0 && y_range ~=0)
    x_plot = y_plot;
elseif(y_range == 0 && x_range ~=0)
    y_plot = x_plot;
end
xlim(x_plot);
ylim(y_plot);

axis equal;

hold off;
