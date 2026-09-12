```matlab
function AF2MecComp()
clc;
close all;


fprintf('============================================================\n');
fprintf('   ANALISIS ESTADISTICO DE CIFRADO DE IMAGENES\n');
fprintf('   Sistemas Dinamicos No Lineales - Mapa Logistico\n');
fprintf('============================================================\n\n');


nombre_archivo = 'Sabrina.jpg';


if ~isfile(nombre_archivo)
    error('No se encontro el archivo "%s" en la carpeta actual.', ...
        nombre_archivo);
end


img_rgb = imread(nombre_archivo);


if ndims(img_rgb) == 2
    img_rgb = cat(3,img_rgb,img_rgb,img_rgb);
end


if ~isa(img_rgb,'uint8')
    img_rgb = im2uint8(img_rgb);
end


img_rgb = imresize(img_rgb,[512 512]);


img_gray = rgb2gray(img_rgb);


fprintf('Imagen utilizada: %s\n',nombre_archivo);
fprintf('Tamano RGB:  %d x %d x 3\n', ...
    size(img_rgb,1),size(img_rgb,2));


fprintf('Tamano Gris: %d x %d\n\n', ...
    size(img_gray,1),size(img_gray,2));




fprintf('============================================================\n');
fprintf('            ANALISIS EN ESCALA DE GRISES\n');
fprintf('============================================================\n');


r_gray  = 3.999;
x0_gray = 0.45321;


img_cifrada_gray = cifrado_caotico( ...
    img_gray,r_gray,x0_gray);


resultados_gray = analizar_estadisticas( ...
    img_gray, ...
    img_cifrada_gray, ...
    'Escala de Grises');




fprintf('\n============================================================\n');
fprintf('                 ANALISIS IMAGEN RGB\n');
fprintf('============================================================\n');


img_cifrada_rgb = zeros(size(img_rgb),'uint8');


img_cifrada_rgb(:,:,1) = cifrado_caotico( ...
    img_rgb(:,:,1), ...
    3.999, ...
    0.12345);


img_cifrada_rgb(:,:,2) = cifrado_caotico( ...
    img_rgb(:,:,2), ...
    3.998, ...
    0.45678);


img_cifrada_rgb(:,:,3) = cifrado_caotico( ...
    img_rgb(:,:,3), ...
    3.997, ...
    0.78912);




[resultados_rgb,resultados_corr_rgb] = ...
    analizar_estadisticas_rgb( ...
    img_rgb, ...
    img_cifrada_rgb);




fprintf('\n\n');
fprintf('============================================================\n');
fprintf('              TABLA - ESCALA DE GRISES\n');
fprintf('============================================================\n');


disp(resultados_gray);


fprintf('\n============================================================\n');
fprintf('                  TABLA - RGB\n');
fprintf('============================================================\n');


disp(resultados_rgb);


fprintf('\n============================================================\n');
fprintf('             CORRELACIONES - RGB\n');
fprintf('============================================================\n');


disp(resultados_corr_rgb);




writetable(resultados_gray,'Resultados_Grises.csv');
writetable(resultados_rgb,'Resultados_RGB.csv');
writetable(resultados_corr_rgb,'Correlaciones_RGB.csv');


fprintf('\nArchivos CSV generados:\n');
fprintf('  - Resultados_Grises.csv\n');
fprintf('  - Resultados_RGB.csv\n');
fprintf('  - Correlaciones_RGB.csv\n');


imwrite(img_cifrada_gray,'Imagen_Cifrada_Grises.png');
imwrite(img_cifrada_rgb,'Imagen_Cifrada_RGB.png');


fprintf('\nImagenes cifradas guardadas:\n');
fprintf('  - Imagen_Cifrada_Grises.png\n');
fprintf('  - Imagen_Cifrada_RGB.png\n');




generar_conclusion( ...
    resultados_gray, ...
    resultados_rgb, ...
    resultados_corr_rgb);


fprintf('\n============================================================\n');
fprintf('                  PROCESO TERMINADO\n');
fprintf('============================================================\n');


end




function img_cifrada = cifrado_caotico(img,r,x0)


[rows,cols] = size(img);


N = rows * cols;


transitorio = 1000;


X = zeros(N + transitorio,1);


X(1) = x0;


for k = 1:(N + transitorio - 1)


    X(k+1) = r * X(k) * (1 - X(k));


end




X = X(transitorio+1 : transitorio+N);


if length(X) ~= N


    error(['La secuencia caotica tiene un numero incorrecto ' ...
           'de elementos.']);


end




U = (2/pi) * asin(sqrt(X));




U = max(U,0);


U = min(U,1-eps);




[~,perm_idx] = sort(U);




img_vector = img(:);




if length(img_vector) ~= N


    error('El numero de pixeles de la imagen no coincide.');


end


if max(perm_idx) > N


    error('La permutacion contiene indices fuera de rango.');


end




img_perm = img_vector(perm_idx);




key_stream = uint8(floor(U * 256));


key_stream = key_stream(:);




cifrado = zeros(N,1,'uint8');




cifrado(1) = bitxor( ...
    img_perm(1), ...
    key_stream(1));


% Resto de pixeles
for k = 2:N


    cifrado(k) = bitxor( ...
        bitxor(img_perm(k),cifrado(k-1)), ...
        key_stream(k));


end




img_cifrada = reshape(cifrado,[rows cols]);


end




function resultados = analizar_estadisticas( ...
    orig,cifrada,titulo)




H_orig = calcular_entropia(orig);


H_cif = calcular_entropia(cifrada);


hist_orig = imhist(orig,256);


hist_cif = imhist(cifrada,256);


CV_orig = calcular_uniformidad(hist_orig);


CV_cif = calcular_uniformidad(hist_cif);


corr_orig = calcular_correlaciones(orig);


corr_cif = calcular_correlaciones(cifrada);




resultados = table( ...
    H_orig, ...
    H_cif, ...
    CV_orig, ...
    CV_cif, ...
    corr_orig.horizontal, ...
    corr_cif.horizontal, ...
    corr_orig.vertical, ...
    corr_cif.vertical, ...
    corr_orig.diagonal, ...
    corr_cif.diagonal, ...
    'VariableNames',{ ...
    'Entropia_Original', ...
    'Entropia_Cifrada', ...
    'CV_Hist_Original', ...
    'CV_Hist_Cifrada', ...
    'Corr_H_Original', ...
    'Corr_H_Cifrada', ...
    'Corr_V_Original', ...
    'Corr_V_Cifrada', ...
    'Corr_D_Original', ...
    'Corr_D_Cifrada'});


fprintf('\nENTROPIA\n');


fprintf('Imagen original : %.6f bits\n',H_orig);


fprintf('Imagen cifrada  : %.6f bits\n',H_cif);


fprintf('\nUNIFORMIDAD DEL HISTOGRAMA\n');


fprintf('CV original : %.6f\n',CV_orig);


fprintf('CV cifrado  : %.6f\n',CV_cif);


fprintf('\nCORRELACION DE PIXELES\n');


fprintf(['Horizontal -> Original: %.6f | ' ...
         'Cifrada: %.6f\n'], ...
    corr_orig.horizontal, ...
    corr_cif.horizontal);


fprintf(['Vertical   -> Original: %.6f | ' ...
         'Cifrada: %.6f\n'], ...
    corr_orig.vertical, ...
    corr_cif.vertical);


fprintf(['Diagonal   -> Original: %.6f | ' ...
         'Cifrada: %.6f\n'], ...
    corr_orig.diagonal, ...
    corr_cif.diagonal);


figure( ...
    'Name',['Analisis - ' titulo], ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Position',[50 50 1400 750]);


subplot(2,3,1);


imshow(orig);


title('Imagen Original','FontWeight','bold');




subplot(2,3,2);


bar(0:255,hist_orig, ...
    'FaceColor',[0.1 0.45 0.80], ...
    'EdgeColor','none');


xlim([0 255]);


grid on;


title('Histograma Original');


xlabel('Nivel de gris');


ylabel('Frecuencia');


subplot(2,3,3);


bar([H_orig,H_cif], ...
    'FaceColor',[0.25 0.65 0.35]);


hold on;


yline(8,'r--','8 bits', ...
    'LineWidth',1.5);


hold off;


ylim([0 8.2]);


grid on;


set(gca,'XTick',[1 2]);


set(gca,'XTickLabel', ...
    {'Original','Cifrada'});


ylabel('Entropia (bits)');


title('Entropia de Shannon');




subplot(2,3,4);


imshow(cifrada);


title('Imagen Cifrada','FontWeight','bold');




subplot(2,3,5);


bar(0:255,hist_cif, ...
    'FaceColor',[0.85 0.25 0.15], ...
    'EdgeColor','none');


xlim([0 255]);


grid on;


title('Histograma Cifrado');


xlabel('Nivel de gris');


ylabel('Frecuencia');




subplot(2,3,6);


datos = [ ...
    corr_orig.horizontal corr_cif.horizontal;
    corr_orig.vertical   corr_cif.vertical;
    corr_orig.diagonal   corr_cif.diagonal];


bar(datos);


hold on;


yline(0,'k-');


hold off;


ylim([-1 1]);


grid on;


set(gca,'XTick',1:3);


set(gca,'XTickLabel', ...
    {'Horizontal','Vertical','Diagonal'});


ylabel('Coeficiente r');


title('Correlacion de Pixeles');


legend('Original','Cifrada', ...
    'Location','best');


nombre_figura = ['Analisis_' ...
    strrep(titulo,' ','_') '.png'];


saveas(gcf,nombre_figura);


end




function [tabla_rgb,tabla_corr] = ...
    analizar_estadisticas_rgb(orig,cifrada)


nombres = {'Red','Green','Blue'};


H_orig = zeros(3,1);
H_cif = zeros(3,1);


CV_orig = zeros(3,1);
CV_cif = zeros(3,1);


H_corr_orig = zeros(3,1);
H_corr_cif = zeros(3,1);


V_corr_orig = zeros(3,1);
V_corr_cif = zeros(3,1);


D_corr_orig = zeros(3,1);
D_corr_cif = zeros(3,1);




for c = 1:3


    canal_orig = orig(:,:,c);


    canal_cif = cifrada(:,:,c);


    H_orig(c) = calcular_entropia(canal_orig);


    H_cif(c) = calcular_entropia(canal_cif);




    h_orig = imhist(canal_orig,256);


    h_cif = imhist(canal_cif,256);


    CV_orig(c) = calcular_uniformidad(h_orig);


    CV_cif(c) = calcular_uniformidad(h_cif);


 
    corr_o = calcular_correlaciones(canal_orig);


    corr_c = calcular_correlaciones(canal_cif);


    H_corr_orig(c) = corr_o.horizontal;
    H_corr_cif(c) = corr_c.horizontal;


    V_corr_orig(c) = corr_o.vertical;
    V_corr_cif(c) = corr_c.vertical;


    D_corr_orig(c) = corr_o.diagonal;
    D_corr_cif(c) = corr_c.diagonal;


end




tabla_rgb = table( ...
    nombres', ...
    H_orig, ...
    H_cif, ...
    CV_orig, ...
    CV_cif, ...
    'VariableNames',{ ...
    'Canal', ...
    'Entropia_Original', ...
    'Entropia_Cifrada', ...
    'CV_Hist_Original', ...
    'CV_Hist_Cifrada'});


tabla_corr = table( ...
    nombres', ...
    H_corr_orig, ...
    H_corr_cif, ...
    V_corr_orig, ...
    V_corr_cif, ...
    D_corr_orig, ...
    D_corr_cif, ...
    'VariableNames',{ ...
    'Canal', ...
    'Horizontal_Original', ...
    'Horizontal_Cifrada', ...
    'Vertical_Original', ...
    'Vertical_Cifrada', ...
    'Diagonal_Original', ...
    'Diagonal_Cifrada'});


fprintf('\nENTROPIA POR CANAL RGB\n');


for c = 1:3


    fprintf('%s -> Original: %.6f bits | Cifrada: %.6f bits\n', ...
        nombres{c}, ...
        H_orig(c), ...
        H_cif(c));


end




fprintf('\nCORRELACIONES RGB\n');


for c = 1:3


    fprintf('\nCanal %s:\n',nombres{c});


    fprintf('  Horizontal -> Original: %.6f | Cifrada: %.6f\n', ...
        H_corr_orig(c),H_corr_cif(c));


    fprintf('  Vertical   -> Original: %.6f | Cifrada: %.6f\n', ...
        V_corr_orig(c),V_corr_cif(c));


    fprintf('  Diagonal   -> Original: %.6f | Cifrada: %.6f\n', ...
        D_corr_orig(c),D_corr_cif(c));


end


figure( ...
    'Name','Analisis Estadistico RGB', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Position',[20 20 1500 850]);


colores = [ ...
    0.85 0.15 0.15;
    0.15 0.65 0.20;
    0.15 0.35 0.85];


subplot(4,4,1);


imshow(orig);


title('Imagen RGB Original', ...
    'FontWeight','bold');




subplot(4,4,9);


imshow(cifrada);


title('Imagen RGB Cifrada', ...
    'FontWeight','bold');




for c = 1:3


    canal_orig = orig(:,:,c);


    canal_cif = cifrada(:,:,c);


    h_orig = imhist(canal_orig,256);


    h_cif = imhist(canal_cif,256);


    subplot(4,4,c+1);


    bar(0:255,h_orig, ...
        'FaceColor',colores(c,:), ...
        'EdgeColor','none');


    xlim([0 255]);


    grid on;


    title(['Canal ' nombres{c} ...
        ' - Original']);


    xlabel('Nivel');


    ylabel('Frecuencia');




    subplot(4,4,c+9);


    bar(0:255,h_cif, ...
        'FaceColor',colores(c,:), ...
        'EdgeColor','none');


    xlim([0 255]);


    grid on;


    title(['Canal ' nombres{c} ...
        ' - Cifrado']);


    xlabel('Nivel');


    ylabel('Frecuencia');


end




subplot(4,4,5);


bar([H_orig,H_cif]);


hold on;


yline(8,'r--','8 bits', ...
    'LineWidth',1.5);


hold off;


ylim([0 8.2]);


grid on;


set(gca,'XTick',1:3);


set(gca,'XTickLabel',nombres);


ylabel('Bits');


title('Entropia por Canal');


legend('Original','Cifrada', ...
    'Location','best');


subplot(4,4,6);


bar([H_corr_orig,H_corr_cif]);


ylim([-1 1]);


grid on;


set(gca,'XTick',1:3);


set(gca,'XTickLabel',nombres);


ylabel('r');


title('Correlacion Horizontal');


legend('Original','Cifrada', ...
    'Location','best');




subplot(4,4,7);


bar([V_corr_orig,V_corr_cif]);


ylim([-1 1]);


grid on;


set(gca,'XTick',1:3);


set(gca,'XTickLabel',nombres);


ylabel('r');


title('Correlacion Vertical');


legend('Original','Cifrada', ...
    'Location','best');




subplot(4,4,8);


bar([D_corr_orig,D_corr_cif]);


ylim([-1 1]);


grid on;


set(gca,'XTick',1:3);


set(gca,'XTickLabel',nombres);


ylabel('r');


title('Correlacion Diagonal');


legend('Original','Cifrada', ...
    'Location','best');




saveas(gcf,'Analisis_RGB.png');


end


function H = calcular_entropia(img)


% Obtener histograma de 256 niveles
histograma = imhist(img,256);


% Convertir a probabilidades
p = histograma / numel(img);


% Eliminar probabilidades cero
p = p(p > 0);


% Entropia de Shannon
H = -sum(p .* log2(p));


end




function CV = calcular_uniformidad(histograma)


h = double(histograma);


media = mean(h);


if media == 0


    CV = Inf;


else


    CV = std(h) / media;


end


end




function resultados = calcular_correlaciones(img)


N_samples = 5000;


[rows,cols] = size(img);




filas = randi(rows,N_samples,1);


columnas = randi(cols-1,N_samples,1);


x = double(img( ...
    sub2ind([rows cols], ...
    filas,columnas)));


y = double(img( ...
    sub2ind([rows cols], ...
    filas,columnas+1)));


R = corrcoef(x,y);


r_horizontal = R(1,2);


filas = randi(rows-1,N_samples,1);


columnas = randi(cols,N_samples,1);


x = double(img( ...
    sub2ind([rows cols], ...
    filas,columnas)));


y = double(img( ...
    sub2ind([rows cols], ...
    filas+1,columnas)));


R = corrcoef(x,y);


r_vertical = R(1,2);




filas = randi(rows-1,N_samples,1);


columnas = randi(cols-1,N_samples,1);


x = double(img( ...
    sub2ind([rows cols], ...
    filas,columnas)));


y = double(img( ...
    sub2ind([rows cols], ...
    filas+1,columnas+1)));


R = corrcoef(x,y);


r_diagonal = R(1,2);


resultados.horizontal = r_horizontal;


resultados.vertical = r_vertical;


resultados.diagonal = r_diagonal;


end




function generar_conclusion( ...
    tabla_gray, ...
    tabla_rgb, ...
    tabla_corr_rgb)


fprintf('\n\n');


fprintf('============================================================\n');
fprintf('                    CONCLUSION\n');
fprintf('============================================================\n\n');


umbral_entropia = 7.90;


umbral_corr = 0.10;


umbral_cv = 0.15;




cumple_entropia_gray = ...
    tabla_gray.Entropia_Cifrada >= umbral_entropia;


cumple_corr_gray = ...
    abs(tabla_gray.Corr_H_Cifrada) <= umbral_corr && ...
    abs(tabla_gray.Corr_V_Cifrada) <= umbral_corr && ...
    abs(tabla_gray.Corr_D_Cifrada) <= umbral_corr;


cumple_uniformidad_gray = ...
    tabla_gray.CV_Hist_Cifrada <= umbral_cv;


cumple_entropia_rgb = ...
    all(tabla_rgb.Entropia_Cifrada >= umbral_entropia);


cumple_uniformidad_rgb = ...
    all(tabla_rgb.CV_Hist_Cifrada <= umbral_cv);


cumple_corr_rgb = ...
    all(abs(tabla_corr_rgb.Horizontal_Cifrada) <= umbral_corr) && ...
    all(abs(tabla_corr_rgb.Vertical_Cifrada) <= umbral_corr) && ...
    all(abs(tabla_corr_rgb.Diagonal_Cifrada) <= umbral_corr);


fprintf('ESCALA DE GRISES\n');


fprintf('  Entropia >= %.2f bits       : %s\n', ...
    umbral_entropia, ...
    estado(cumple_entropia_gray));


fprintf('  Histograma uniforme         : %s\n', ...
    estado(cumple_uniformidad_gray));


fprintf('  Correlacion |r| <= %.2f     : %s\n', ...
    umbral_corr, ...
    estado(cumple_corr_gray));


fprintf('\nIMAGEN RGB\n');


fprintf('  Entropia >= %.2f bits       : %s\n', ...
    umbral_entropia, ...
    estado(cumple_entropia_rgb));


fprintf('  Histogramas uniformes       : %s\n', ...
    estado(cumple_uniformidad_rgb));


fprintf('  Correlacion |r| <= %.2f     : %s\n', ...
    umbral_corr, ...
    estado(cumple_corr_rgb));


fprintf('\n------------------------------------------------------------\n');


if cumple_entropia_gray && ...
        cumple_corr_gray && ...
        cumple_entropia_rgb && ...
        cumple_corr_rgb


    fprintf('RESULTADO GENERAL: SATISFACTORIO\n\n');


    fprintf(['El proceso de cifrado basado en el mapa logistico ' ...
        'caotico presenta un comportamiento adecuado desde ' ...
        'el punto de vista estadistico. La entropia de las ' ...
        'imagenes cifradas se aproxima al valor ideal de 8 ' ...
        'bits, indicando un incremento de la incertidumbre ' ...
        'en la distribucion de los niveles de intensidad. ' ...
        'Ademas, la correlacion entre pixeles adyacentes se ' ...
        'reduce considerablemente respecto a las imagenes ' ...
        'originales y se aproxima a cero en las direcciones ' ...
        'horizontal, vertical y diagonal. Esto indica que el ' ...
        'cifrado disminuye la dependencia espacial presente ' ...
        'en la imagen original.\n']);


else


    fprintf('RESULTADO GENERAL: PARCIALMENTE SATISFACTORIO\n\n');


    fprintf(['Los resultados muestran que el proceso de cifrado ' ...
        'reduce la dependencia estadistica entre pixeles y ' ...
        'aumenta la incertidumbre de la imagen. Sin embargo, ' ...
        'uno o mas criterios establecidos no alcanzaron ' ...
        'completamente los valores ideales. En particular, ' ...
        'puede ser necesario mejorar la etapa de difusion ' ...
        'para obtener una distribucion mas uniforme de los ' ...
        'niveles de intensidad y una entropia mas cercana a ' ...
        '8 bits.\n']);


end


fprintf('------------------------------------------------------------\n');


end


%% ################################################################


function texto = estado(condicion)


if condicion


    texto = 'CUMPLE';


else


    texto = 'NO CUMPLE';


end


end


```
