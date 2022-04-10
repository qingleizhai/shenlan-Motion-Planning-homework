clear
clc
syms T
syms p_x_f p_y_f p_z_f v_x_0 v_y_0 v_z_0 p_x_0 p_y_0 p_z_0 v_x_f v_y_f v_z_f

%solve for coeffs:
 A = [
        -12/(T^3)         0         0  6/(T^2)       0            0; 
                0 -12/(T^3)         0        0  6/(T^2)           0; 
                0         0 -12/(T^3)        0        0     6/(T^2); 
          6/(T^2)         0         0     -2/T        0           0; 
                0   6/(T^2)         0        0     -2/T           0; 
                0         0    6/(T^2)       0        0        -2/T
      ];
 B = [
        p_x_f - v_x_0 * T - p_x_0; 
        p_y_f - v_y_0 * T - p_y_0;
        p_z_f - v_z_0 * T - p_z_0;
                    v_x_f - v_x_0;
                    v_y_f - v_y_0;
                    v_z_f - v_z_0
      ];
c = simplify(A * B);

%cost function
M = [
        (T^3)/3       0       0 (T^2)/2        0        0;
              0 (T^3)/3       0       0  (T^2)/2        0;
              0       0 (T^3)/3       0        0  (T^2)/2;
        (T^2)/2       0       0       T        0        0;
              0 (T^2)/2       0       0        T        0;
              0       0 (T^2)/2       0        0        T
    ];
J = collect(expand(simplify(transpose(c) * M * c) + T),T);

%derivative
d = T^3;
n = d * J;
der = collect(expand(simplify(diff(n, T) * d - n * diff(d, T))),T);

[c,t] = coeffs(der,T)


