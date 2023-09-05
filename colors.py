import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Carregar a tabela de sensibilidade da CIE (p.csv)
cie_table = pd.read_csv('p.csv', skiprows=1, delimiter=',', names=['l', 'x', 'y', 'z'])

# Definir a temperatura do sol em Kelvin
temperature = 5778  # Temperatura do sol

# Definir a faixa de comprimentos de onda relevantes (por exemplo, de 380 nm a 780 nm)
wavelengths = np.linspace(380, 780, 401)  # 401 pontos entre 380 e 780 nm
delta_lambda = wavelengths[1] - wavelengths[0]

# Cálculo da radiação espectral do sol com dispersão de Rayleigh
h = 6.626e-34  # Constante de Planck em J*s
c = 3e8  # Velocidade da luz em m/s
k = 1.38e-23  # Constante de Boltzmann em J/K

radiance = []
for wavelength in wavelengths:
    # Comprimento de onda em metros
    lambda_m = wavelength * 1e-9
    # Cálculo da radiância espectral do sol com dispersão de Rayleigh
    if temperature > 0:
        radiance_value = (2 * h * c ** 2 / (lambda_m ** 5)) / (np.exp(h * c / (lambda_m * k * temperature)) - 1)
        radiance.append(radiance_value if not np.isnan(radiance_value) else 0)
    else:
        radiance.append(0)

radiance = np.array(radiance)

# Interpolação dos valores da tabela de sensibilidade da CIE para os comprimentos de onda correspondentes
cie_x = np.interp(wavelengths, cie_table['l'], cie_table['x'])
cie_y = np.interp(wavelengths, cie_table['l'], cie_table['y'])
cie_z = np.interp(wavelengths, cie_table['l'], cie_table['z'])

# Calcular os valores ponderados da sensibilidade da CIE
x_value = np.nansum(radiance * cie_x) * delta_lambda
y_value = np.nansum(radiance * cie_y) * delta_lambda
z_value = np.nansum(radiance * cie_z) * delta_lambda

# Converter para coordenadas XYZ
xyz_sum = x_value + y_value + z_value
x_normalized = x_value / xyz_sum
y_normalized = y_value / xyz_sum
z_normalized = z_value / xyz_sum

# Converter para coordenadas RGB (D65 Illuminant)
rgb_r = x_normalized * 3.2406 + y_normalized * -1.5372 + z_normalized * -0.4986
rgb_g = x_normalized * -0.9689 + y_normalized * 1.8758 + z_normalized * 0.0415
rgb_b = x_normalized * 0.0557 + y_normalized * -0.2040 + z_normalized * 1.0570

# Ajustar os valores de RGB para o intervalo [0, 1]
max_value = max(rgb_r, rgb_g, rgb_b)
rgb_r /= max_value if max_value > 0 else 1
rgb_g /= max_value if max_value > 0 else 1
rgb_b /= max_value if max_value > 0 else 1

# Imprimir os valores RGB
print(f"Valores RGB ajustados: (R = {rgb_r:.4f}, G = {rgb_g:.4f}, B = {rgb_b:.4f}")

# Plotar a cor percebida
color = (rgb_r, rgb_g, rgb_b)
plt.imshow([[color]])
plt.axis('off')
plt.show()

