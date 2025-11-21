clc; clear; close all;

textura=imread("C:\Users\Simone\Downloads\textura.jpg");

% volcan
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
r = sqrt(x.^2 + y.^2);

z = exp(-0.1*r.^2) * 5 - 2*exp(-r.^2);
z = z + 0.07*randn(size(z));

surf(x*1000, y*1000, z*1000,'CData',textura,'FaceColor','texturemap','EdgeColor','none');
axis equal
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Altura (m)')
title('Modelo 3D de un volcán con erupciones')
lighting phong
camlight headlight
hold on

% constantes
g = 9.81;
h0 = 4100;   % altura del volcán
k = 0.003;   % arrastre


%viento
v_wx = 10+rand*4;    % viento base en X
v_wy = -5+rand*4;    % viento base en Y

gamma_wind = 0.0003;   % incremento del viento con la altura


%magnus
omega = 30;         % velocidad angular
alphaM = 0.0005;     % intensidad Magnus


% colores
colores = [
    251 116 168;
    252 172 57;
    251 149 1;
    218 43 66;
    245 77 233
] / 255;

N = 15;
v = [50, 75, 100, 125, 150, 175, 200];


% erupciones
for e = 1:5
    trayectorias = struct();
    max_t = 0;

    for i = 1:N
        % Velocidades iniciales
        v0 = v(randi(length(v)));
        theta = rand * 25 + 20;
        phi = rand * 360;

        vx = v0 * cosd(theta) * cosd(phi);
        vy = v0 * cosd(theta) * sind(phi);
        vz = v0 * sind(theta);

        % Tiempo de vuelo (ajustado por arrastre)
        coef = [0.5*g, -vz, -h0];
        t_sol = roots(coef);
        t_vuelo = max(t_sol);

        f = @(t) h0 + (vz + g/k)/k * (1 - exp(-k*t)) - g*t/k;
        t_vuelo_r = fzero(f, t_vuelo);

        t = linspace(0, t_vuelo_r, 200);

        % trayectoria
        X = zeros(size(t));
        Y = zeros(size(t));
        Z = zeros(size(t));

        for k2 = 1:length(t)

            % Movimiento vertical
            Z(k2) = h0 + (vz + g/k)/k * (1 - exp(-k*t(k2))) - g*t(k2)/k;

            % Viento dependiente de altura
            vwx_z = v_wx * (1 + gamma_wind * Z(k2));
            vwy_z = v_wy * (1 + gamma_wind * Z(k2));

            % Velocidad horizontal 
            vx_eff = vx + vwx_z;
            vy_eff = vy + vwy_z;

            % Movimiento horizontal con arrastre
            X(k2) = (vx_eff/k) * (1 - exp(-k*t(k2)));
            Y(k2) = (vy_eff/k) * (1 - exp(-k*t(k2)));

            % Efecto Magnus basado en velocidad
            X(k2) = X(k2) + alphaM * omega * vy_eff * (1 - exp(-k*t(k2)));
            Y(k2) = Y(k2) - alphaM * omega * vx_eff * (1 - exp(-k*t(k2)));

        end

        % dibujar
        color = colores(mod(i-1, size(colores,1)) + 1, :);

        trayectorias(i).X = X;
        trayectorias(i).Y = Y;
        trayectorias(i).Z = Z;
        trayectorias(i).color = color;

        trayectorias(i).marker = plot3(X(1), Y(1), Z(1), 'o', ...
            'Color', color, 'MarkerFaceColor', color, ...
            'MarkerEdgeColor', [0,0,0], 'LineWidth', 1);

        trayectorias(i).line = plot3(nan, nan, nan, '--', ...
            'Color', color, 'LineWidth', 0.7);

        max_t = max(max_t, length(t));
    end


    % animacion
    for j = 1:max_t
        for i = 1:N
            X = trayectorias(i).X;
            Y = trayectorias(i).Y;
            Z = trayectorias(i).Z;

            if j <= length(X)
                set(trayectorias(i).marker, 'XData', X(j), 'YData', Y(j), 'ZData', Z(j));
                set(trayectorias(i).line, 'XData', X(1:j), 'YData', Y(1:j), 'ZData', Z(1:j));
            end
        end

        drawnow limitrate
        pause(0.01)
    end

    pause(1);
end

view(45, 30);