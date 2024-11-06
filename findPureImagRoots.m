function [Gain, w_j] = findPureImagRoots(Ds)

%   ����:
%       Ds - ����ʽϵ���������Ӹߴε��ʹ����С�����:
%            Ds = [1, 3, 12, K_sym - 16, K_sym]
%
%   ���:
%       params      - �ṹ�壬�������Ų����Ľ�
%       frequencies - ��Ӧ��Ƶ�� w �Ľ⣨��������Ƶ�ʣ�
%
%   
%     

    % ��ȡ����ʽ�ķ��ű���
    vars = symvar(Ds);

    % ��һ������ʽ��ʹ��ߴ���ϵ��Ϊ 1
    if Ds(1) ~= 1
        Ds = Ds / Ds(1);
        fprintf('����ʽ�ѹ�һ������ߴ���ϵ������Ϊ 1��\n');
    end

    syms w real
    syms s
    syms K_x real
    % �����Ϊ������ s = jw

    % ������ʽ���� s = jw
    Ds_sub = poly2sym(Ds, s);
    Ds_sub = subs(Ds_sub, s,1i*w);
    Ds_sub = subs(Ds_sub,vars,K_x);
    % չ������ʽ
    Ds_expanded = expand(Ds_sub);

    % ����ʵ�����鲿
    Ds_real = real(Ds_expanded);
    Ds_imag = imag(Ds_expanded);

    % ���������飺ʵ�� = 0 ���鲿 = 0
    eq1 = Ds_real == 0;
    eq2 = Ds_imag == 0;

    % �ϲ�����
    eqs = [eq1, eq2];

    % ��ⷽ����
    [w_j, Gain] = solve(eqs, [w, K_x], 'Real', true);

    % ����ת��Ϊ��ֵ��ʽ
    w_j = double(w_j);
    
    % ����Ƿ��н�
    if isempty(w_j)
        warning('δ�ҵ�ʹ����ʽ��Ϊ�������Ĳ���ֵ��');
        Gain = [];
        w_j = [];
        return;
    end
    Gain = double(Gain);
    maski = (w_j == 0 & Gain == 0) | Gain < 0;
    w_j(maski)  = [];
    Gain(maski) = [];
end





