function [Gain, w_j] = findPureImagRoots(Ds)

%   输入:
%       Ds - 多项式系数向量，从高次到低次排列。例如:
%            Ds = [1, 3, 12, K_sym - 16, K_sym]
%
%   输出:
%       params      - 结构体，包含符号参数的解
%       frequencies - 对应的频率 w 的解（包括正负频率）
%
%   
%     

    % 获取多项式的符号变量
    vars = symvar(Ds);

    % 归一化多项式，使最高次项系数为 1
    if Ds(1) ~= 1
        Ds = Ds / Ds(1);
        fprintf('多项式已归一化，最高次项系数设置为 1。\n');
    end

    syms w real
    syms s
    syms K_x real
    % 假设根为纯虚数 s = jw

    % 将多项式代入 s = jw
    Ds_sub = poly2sym(Ds, s);
    Ds_sub = subs(Ds_sub, s,1i*w);
    Ds_sub = subs(Ds_sub,vars,K_x);
    % 展开多项式
    Ds_expanded = expand(Ds_sub);

    % 分离实部和虚部
    Ds_real = real(Ds_expanded);
    Ds_imag = imag(Ds_expanded);

    % 建立方程组：实部 = 0 和虚部 = 0
    eq1 = Ds_real == 0;
    eq2 = Ds_imag == 0;

    % 合并方程
    eqs = [eq1, eq2];

    % 求解方程组
    [w_j, Gain] = solve(eqs, [w, K_x], 'Real', true);

    % 将解转换为数值形式
    w_j = double(w_j);
    
    % 检查是否有解
    if isempty(w_j)
        warning('未找到使多项式根为纯虚数的参数值。');
        Gain = [];
        w_j = [];
        return;
    end
    Gain = double(Gain);
    maski = (w_j == 0 & Gain == 0) | Gain < 0;
    w_j(maski)  = [];
    Gain(maski) = [];
end





